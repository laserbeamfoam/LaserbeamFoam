#!/usr/bin/env python3
"""GUI case spawner for laserbeamFoam tutorial templates.

The tool is built around three ideas:
1. Walk a tutorial case and extract editable entries from OpenFOAM-style inputs.
2. Define sweep axes against any extracted entry and preview the cartesian array.
3. Materialize the generated variants into a clean output tree without copying run outputs.

It also reuses the shared tutorials/plot_openfoam_log.py parser to summarize existing
laserbeamFoam log files and treats dense numeric assets like beamShape.inp as previewable
matrix fields, with optional CuPy acceleration for preview normalization.

Typical usage:
	conda run -n ds python spawn.py /path/to/tutorials/beam_profiles/laserbeamFoam/opa_cbc_movingFrame_ss316L
"""

from __future__ import annotations

import importlib.util
import itertools
import json
import math
import os
import re
import shutil
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np

try:
	import cupy as cp

	HAVE_CUPY = True
except Exception:
	cp = None
	HAVE_CUPY = False

try:
	from OpenGL import GL, GLU

	HAVE_PYOPENGL = True
except Exception:
	GL = None
	GLU = None
	HAVE_PYOPENGL = False

from PySide6.QtCore import QPoint, Qt
from PySide6.QtGui import QAction, QColor, QImage, QPainter, QPixmap, QSurfaceFormat
from PySide6.QtOpenGLWidgets import QOpenGLWidget
from PySide6.QtWidgets import (
	QApplication,
	QAbstractItemView,
	QComboBox,
	QFileDialog,
	QFormLayout,
	QGroupBox,
	QHBoxLayout,
	QHeaderView,
	QLabel,
	QLineEdit,
	QMainWindow,
	QMessageBox,
	QPlainTextEdit,
	QProgressDialog,
	QPushButton,
	QSizePolicy,
	QSplitter,
	QStatusBar,
	QTableWidget,
	QTableWidgetItem,
	QToolBar,
	QTreeWidget,
	QTreeWidgetItem,
	QVBoxLayout,
	QWidget,
)


THIS_FILE = Path(__file__).resolve()
SHARED_PLOTTER = THIS_FILE.parents[2] / "plot_openfoam_log.py"

TIME_DIR_RE = re.compile(r"^\d+(?:\.\d+)?$")
PROCESSOR_DIR_RE = re.compile(r"^processor\d+$")
FLOAT_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?$")
NUMERIC_LINE_RE = re.compile(r"^\s*[-+0-9.eE\s]+\s*$")
TRAILING_OUTPUT_SUFFIXES = {".foam", ".png", ".jpg", ".jpeg", ".gif", ".svg", ".vtk"}
SKIP_DIR_NAMES = {"postProcessing", "VTK", "__pycache__", ".git"}
PUNCTUATION = {"{", "}", "(", ")", ";", "[", "]"}


@dataclass(slots=True)
class CaseEntry:
	relative_path: Path
	path_segments: tuple[str, ...]
	value_start: int
	value_end: int
	value_text: str

	@property
	def path_label(self) -> str:
		return "/".join(self.path_segments)

	@property
	def identifier(self) -> str:
		return f"{self.relative_path.as_posix()}::{self.path_label}"

	@property
	def short_label(self) -> str:
		if not self.path_segments:
			return self.relative_path.name
		if len(self.path_segments) == 1:
			return self.path_segments[0]
		return f"{self.path_segments[-2]}/{self.path_segments[-1]}"

	@property
	def display_value(self) -> str:
		return compact_text(self.value_text)


@dataclass(slots=True)
class ParsedFile:
	path: Path
	relative_path: Path
	kind: str
	entries: list[CaseEntry] = field(default_factory=list)
	text: str | None = None
	parse_error: str | None = None


@dataclass(slots=True)
class SweepDefinition:
	entry_id: str
	label: str
	mode: str
	spec: str
	enabled: bool = True


@dataclass(slots=True)
class VariantDescriptor:
	index: int
	case_name: str
	assignments: list[tuple[SweepDefinition, str]]


@dataclass(slots=True)
class Token:
	kind: str
	text: str
	start: int
	end: int


class ParseError(RuntimeError):
	pass


class SweepSpecError(RuntimeError):
	pass


def compact_text(text: str, limit: int = 120) -> str:
	single_line = " ".join(text.split())
	if len(single_line) <= limit:
		return single_line
	return single_line[: limit - 3] + "..."


def is_case_directory(path: Path) -> bool:
	return (path / "system").is_dir() and (path / "constant").is_dir()


def default_case_root() -> Path:
	cwd = Path.cwd().resolve()
	if is_case_directory(cwd):
		return cwd

	parent = THIS_FILE.parent
	if is_case_directory(parent):
		return parent

	for child in sorted(parent.iterdir()):
		if child.is_dir() and is_case_directory(child):
			return child

	return parent


def sanitize_token(text: str, limit: int = 48) -> str:
	token = re.sub(r"[^A-Za-z0-9._+-]+", "-", text.strip())
	token = token.strip("-._")
	if not token:
		token = "value"
	if len(token) > limit:
		token = token[:limit].rstrip("-._")
	return token or "value"


def try_parse_float(text: str) -> float | None:
	stripped = text.strip()
	if not FLOAT_RE.match(stripped):
		return None
	try:
		return float(stripped)
	except ValueError:
		return None


def format_numeric(value: float, template: str | None = None) -> str:
	if template:
		stripped = template.strip()
		if stripped and re.fullmatch(r"[-+]?\d+", stripped):
			rounded = int(round(value))
			return str(rounded)
		if "e" in stripped.lower():
			precision = 8
			if "." in stripped:
				mantissa = stripped.lower().split("e", maxsplit=1)[0]
				decimals = mantissa.split(".", maxsplit=1)[1]
				precision = max(1, len(decimals))
			return f"{value:.{precision}e}"
		if "." in stripped:
			decimals = len(stripped.split(".", maxsplit=1)[1])
			if decimals > 0:
				rendered = f"{value:.{decimals}f}"
				return rendered.rstrip("0").rstrip(".") if decimals > 1 else rendered
	rendered = f"{value:.12g}"
	return rendered


def looks_like_openfoam_text(relative_path: Path, text: str) -> bool:
	if "FoamFile" in text:
		return True

	top = relative_path.parts[0] if relative_path.parts else ""
	if top in {"system", "constant", "0", "initial"}:
		return ";" in text or "{" in text

	return False


def looks_like_numeric_matrix(text: str) -> bool:
	lines = [line for line in text.splitlines() if line.strip()]
	if len(lines) < 2:
		return False

	sample = lines[: min(24, len(lines))]
	widths: list[int] = []
	for line in sample:
		if not NUMERIC_LINE_RE.match(line):
			return False
		row = np.fromstring(line, sep=" ")
		if row.size < 4:
			return False
		widths.append(int(row.size))

	return bool(widths) and max(widths) - min(widths) <= 2


def tokenize_text(text: str) -> list[Token]:
	tokens: list[Token] = []
	index = 0
	length = len(text)

	while index < length:
		char = text[index]

		if char.isspace():
			index += 1
			continue

		if text.startswith("//", index):
			newline = text.find("\n", index)
			index = length if newline == -1 else newline + 1
			continue

		if text.startswith("/*", index):
			close = text.find("*/", index + 2)
			index = length if close == -1 else close + 2
			continue

		if char in {'"', "'"}:
			start = index
			quote = char
			index += 1
			while index < length:
				if text[index] == "\\":
					index += 2
					continue
				if text[index] == quote:
					index += 1
					break
				index += 1
			tokens.append(Token("string", text[start:index], start, index))
			continue

		if char in PUNCTUATION:
			tokens.append(Token("symbol", char, index, index + 1))
			index += 1
			continue

		start = index
		if char == "#":
			index += 1
			while index < length and not text[index].isspace() and text[index] not in PUNCTUATION and text[index] not in {'"', "'"}:
				index += 1
			tokens.append(Token("directive", text[start:index], start, index))
			continue

		while index < length:
			if text[index].isspace() or text[index] in PUNCTUATION:
				break
			if text.startswith("//", index) or text.startswith("/*", index):
				break
			index += 1
		tokens.append(Token("word", text[start:index], start, index))

	return tokens


class FoamEntryExtractor:
	def __init__(self, text: str, relative_path: Path):
		self.text = text
		self.relative_path = relative_path
		self.tokens = tokenize_text(text)
		self.index = 0
		self.entries: list[CaseEntry] = []

	def extract(self) -> list[CaseEntry]:
		self._parse_sequence((), None)
		return self.entries

	def _current(self) -> Token | None:
		if self.index >= len(self.tokens):
			return None
		return self.tokens[self.index]

	def _peek_text(self, offset: int = 0) -> str | None:
		position = self.index + offset
		if position >= len(self.tokens):
			return None
		return self.tokens[position].text

	def _expect(self, symbol: str) -> None:
		token = self._current()
		if token is None or token.text != symbol:
			raise ParseError(f"Expected '{symbol}' in {self.relative_path.as_posix()}")

	def _token_value(self, token: Token) -> str:
		if token.kind == "string":
			return token.text[1:-1]
		return token.text

	def _add_entry(self, path_segments: tuple[str, ...], start: int, end: int) -> None:
		value_text = self.text[start:end]
		if not value_text.strip():
			return
		self.entries.append(
			CaseEntry(
				relative_path=self.relative_path,
				path_segments=path_segments,
				value_start=start,
				value_end=end,
				value_text=value_text,
			)
		)

	def _parse_sequence(self, context: tuple[str, ...], end_symbol: str | None) -> None:
		while True:
			token = self._current()
			if token is None:
				return
			if end_symbol is not None and token.text == end_symbol:
				return
			if token.text == ";":
				self.index += 1
				continue
			if token.kind in {"word", "string", "directive"}:
				self._parse_statement(context)
				continue
			if token.text in {"(", "{"}:
				closer = ")" if token.text == "(" else "}"
				self._skip_group(token.text, closer)
				continue
			self.index += 1

	def _parse_statement(self, context: tuple[str, ...]) -> None:
		token = self._current()
		if token is None:
			return

		self.index += 1
		if token.kind == "directive" or token.text.startswith("#"):
			self._skip_directive()
			return

		key = self._token_value(token)
		current = self._current()
		if current is None:
			return

		if current.text == "{":
			self._parse_block(context + (key,))
			return

		if current.text == "(":
			has_children, start, end = self._parse_list(context + (key,))
			if not has_children:
				self._add_entry(context + (key,), start, end)
			if self._peek_text() == ";":
				self.index += 1
			return

		start_end = self._consume_value_span()
		if start_end is not None:
			start, end = start_end
			self._add_entry(context + (key,), start, end)
		if self._peek_text() == ";":
			self.index += 1

	def _parse_block(self, context: tuple[str, ...]) -> None:
		self._expect("{")
		self.index += 1
		self._parse_sequence(context, "}")
		self._expect("}")
		self.index += 1
		if self._peek_text() == ";":
			self.index += 1

	def _parse_list(self, context: tuple[str, ...]) -> tuple[bool, int, int]:
		self._expect("(")
		start = self._current().start
		self.index += 1
		has_children = False

		while True:
			token = self._current()
			if token is None:
				raise ParseError(f"Unterminated list in {self.relative_path.as_posix()}")
			if token.text == ")":
				end = token.end
				self.index += 1
				return has_children, start, end
			if token.kind in {"word", "string"} and self._peek_text(1) == "{":
				has_children = True
				name = self._token_value(token)
				self.index += 1
				self._parse_block(context + (name,))
				continue
			if token.kind == "directive":
				self._skip_directive()
				continue
			if token.text == "(":
				self._skip_group("(", ")")
				continue
			if token.text == "{":
				self._skip_group("{", "}")
				continue
			self.index += 1

	def _consume_value_span(self) -> tuple[int, int] | None:
		token = self._current()
		if token is None:
			return None

		start = token.start
		end = token.end
		paren_depth = 0
		brace_depth = 0
		bracket_depth = 0

		while True:
			token = self._current()
			if token is None:
				break
			if token.text == ";" and paren_depth == 0 and brace_depth == 0 and bracket_depth == 0:
				break
			if token.text == "(" :
				paren_depth += 1
			elif token.text == ")":
				if paren_depth == 0:
					break
				paren_depth -= 1
			elif token.text == "{":
				brace_depth += 1
			elif token.text == "}":
				if brace_depth == 0:
					break
				brace_depth -= 1
			elif token.text == "[":
				bracket_depth += 1
			elif token.text == "]":
				bracket_depth = max(0, bracket_depth - 1)
			end = token.end
			self.index += 1

		return start, end

	def _skip_group(self, opener: str, closer: str) -> None:
		self._expect(opener)
		depth = 0
		while True:
			token = self._current()
			if token is None:
				raise ParseError(f"Unterminated group in {self.relative_path.as_posix()}")
			if token.text == opener:
				depth += 1
			elif token.text == closer:
				depth -= 1
			self.index += 1
			if depth == 0:
				return

	def _skip_directive(self) -> None:
		while True:
			token = self._current()
			if token is None:
				return
			if token.text in {";", "{", "}"}:
				break
			self.index += 1
		if self._peek_text() == ";":
			self.index += 1


def parse_foam_entries(text: str, relative_path: Path) -> list[CaseEntry]:
	extractor = FoamEntryExtractor(text, relative_path)
	return extractor.extract()


def should_skip_directory(relative_path: Path) -> bool:
	if not relative_path.parts:
		return False
	name = relative_path.name
	if name in SKIP_DIR_NAMES:
		return True
	if PROCESSOR_DIR_RE.match(name):
		return True
	if TIME_DIR_RE.match(name) and name != "0":
		return True
	return False


def should_skip_generated_path(relative_path: Path) -> bool:
	name = relative_path.name
	if should_skip_directory(relative_path):
		return True
	if name.startswith("log."):
		return True
	if name in {"plotLog.png", "spawn_manifest.json", "spawn_variant.json"}:
		return True
	if relative_path.suffix.lower() in TRAILING_OUTPUT_SUFFIXES:
		return True
	return False


def scan_case(root: Path) -> dict[str, ParsedFile]:
	parsed_files: dict[str, ParsedFile] = {}

	for dirpath, dirnames, filenames in os.walk(root):
		current_dir = Path(dirpath)
		try:
			rel_dir = current_dir.relative_to(root)
		except ValueError:
			continue

		dirnames[:] = [name for name in dirnames if not should_skip_directory(rel_dir / name)]

		for filename in sorted(filenames):
			relative_path = rel_dir / filename if rel_dir != Path(".") else Path(filename)
			path = root / relative_path

			if path.suffix.lower() in {".zip", ".pyc"}:
				continue

			if filename.startswith("log."):
				parsed_files[relative_path.as_posix()] = ParsedFile(
					path=path,
					relative_path=relative_path,
					kind="log",
				)
				continue

			try:
				text = path.read_text(encoding="utf-8", errors="ignore")
			except OSError as exc:
				parsed_files[relative_path.as_posix()] = ParsedFile(
					path=path,
					relative_path=relative_path,
					kind="unsupported",
					parse_error=str(exc),
				)
				continue

			if looks_like_numeric_matrix(text):
				parsed_files[relative_path.as_posix()] = ParsedFile(
					path=path,
					relative_path=relative_path,
					kind="numeric",
				)
				continue

			if looks_like_openfoam_text(relative_path, text):
				try:
					entries = parse_foam_entries(text, relative_path)
					parsed_files[relative_path.as_posix()] = ParsedFile(
						path=path,
						relative_path=relative_path,
						kind="foam",
						entries=entries,
						text=text,
					)
					continue
				except ParseError as exc:
					parsed_files[relative_path.as_posix()] = ParsedFile(
						path=path,
						relative_path=relative_path,
						kind="text",
						text=text,
						parse_error=str(exc),
					)
					continue

			parsed_files[relative_path.as_posix()] = ParsedFile(
				path=path,
				relative_path=relative_path,
				kind="text",
				text=text,
			)

	return parsed_files


@lru_cache(maxsize=8)
def load_numeric_matrix(path_text: str) -> np.ndarray:
	path = Path(path_text)
	text = path.read_text(encoding="utf-8", errors="ignore")
	lines = [line for line in text.splitlines() if line.strip()]
	rows = [np.fromstring(line, sep=" ", dtype=np.float64) for line in lines]
	width = max((row.size for row in rows), default=0)
	matrix = np.full((len(rows), width), np.nan, dtype=np.float64)
	for index, row in enumerate(rows):
		matrix[index, : row.size] = row
	return matrix


def downsample_matrix(matrix: np.ndarray, max_dim: int = 640) -> np.ndarray:
	row_step = max(1, math.ceil(matrix.shape[0] / max_dim))
	col_step = max(1, math.ceil(matrix.shape[1] / max_dim))
	return matrix[::row_step, ::col_step]


def normalize_matrix_for_preview(matrix: np.ndarray) -> np.ndarray:
	sampled = downsample_matrix(matrix)

	if HAVE_CUPY and sampled.size > 262_144:
		gpu_matrix = cp.asarray(np.nan_to_num(sampled, nan=0.0, posinf=0.0, neginf=0.0))
		min_value = float(cp.min(gpu_matrix).get())
		max_value = float(cp.max(gpu_matrix).get())
		if max_value - min_value < 1.0e-12:
			normalized = cp.zeros_like(gpu_matrix)
		else:
			normalized = cp.clip((gpu_matrix - min_value) / (max_value - min_value), 0.0, 1.0)
		return cp.asnumpy(normalized)

	finite = np.nan_to_num(sampled, nan=0.0, posinf=0.0, neginf=0.0)
	min_value = float(np.min(finite))
	max_value = float(np.max(finite))
	if max_value - min_value < 1.0e-12:
		return np.zeros_like(finite)
	return np.clip((finite - min_value) / (max_value - min_value), 0.0, 1.0)


def matrix_to_qimage(matrix: np.ndarray) -> QImage:
	normalized = normalize_matrix_for_preview(matrix)
	rgb = np.stack(
		[
			normalized,
			np.sqrt(np.clip(normalized, 0.0, 1.0)),
			1.0 - normalized,
		],
		axis=-1,
	)
	rgb8 = np.ascontiguousarray((rgb * 255.0).astype(np.uint8))
	height, width, _ = rgb8.shape
	image = QImage(rgb8.data, width, height, width * 3, QImage.Format.Format_RGB888)
	return image.copy()


def describe_matrix(matrix: np.ndarray) -> str:
	finite = matrix[np.isfinite(matrix)]
	if finite.size == 0:
		return f"shape={matrix.shape[0]} x {matrix.shape[1]} | no finite values"
	source = "CuPy" if HAVE_CUPY and matrix.size > 262_144 else "NumPy"
	return (
		f"shape={matrix.shape[0]} x {matrix.shape[1]} | min={float(np.min(finite)):.6g} | "
		f"max={float(np.max(finite)):.6g} | mean={float(np.mean(finite)):.6g} | preview={source}"
	)


@lru_cache(maxsize=1)
def load_shared_plotter() -> Any:
	if not SHARED_PLOTTER.is_file():
		raise FileNotFoundError(f"Shared log plotter not found: {SHARED_PLOTTER}")

	spec = importlib.util.spec_from_file_location("laserbeamfoam_plot_openfoam_log", SHARED_PLOTTER)
	if spec is None or spec.loader is None:
		raise RuntimeError(f"Could not load module spec for {SHARED_PLOTTER}")

	module = importlib.util.module_from_spec(spec)
	sys.modules[spec.name] = module
	spec.loader.exec_module(module)
	return module


def summarize_log(path: Path) -> str:
	plotter = load_shared_plotter()
	steps = plotter.parse_log(path)
	if not steps:
		return f"No time-step data found in {path.name}."

	latest = steps[-1]
	fields = plotter.discovered_fields(steps)
	duration = latest.execution_time if latest.execution_time is not None else float("nan")
	lines = [
		f"Log file: {path}",
		f"Time steps parsed: {len(steps)}",
		f"Latest simulation time: {latest.time:.6e}",
	]

	if math.isfinite(duration):
		lines.append(f"Latest execution time: {duration:.6f} s")

	if latest.delta_t is not None:
		lines.append(f"Latest deltaT: {latest.delta_t:.6e}")

	if latest.courant_max is not None:
		lines.append(f"Latest Co max: {latest.courant_max:.6g}")

	if latest.laser_power is not None:
		lines.append(f"Latest laser power: {latest.laser_power:.6g}")

	if latest.total_q_deposited is not None:
		lines.append(f"Latest deposited energy: {latest.total_q_deposited:.6g}")

	if fields:
		preview_fields = ", ".join(fields[:12])
		lines.append(f"Discovered fields: {preview_fields}")
		if len(fields) > 12:
			lines.append(f"Additional fields: {len(fields) - 12}")

	return "\n".join(lines)


def split_list_spec(spec: str) -> list[str]:
	items: list[str] = []
	for raw in re.split(r"[\n,]", spec):
		value = raw.strip()
		if value:
			items.append(value)
	return items


def evaluate_sweep(definition: SweepDefinition, template_value: str) -> list[str]:
	if not definition.enabled:
		return []

	spec = definition.spec.strip()
	if not spec:
		raise SweepSpecError(f"Sweep '{definition.label}' has an empty specification.")

	if definition.mode == "list":
		values = split_list_spec(spec)
		if not values:
			raise SweepSpecError(f"Sweep '{definition.label}' does not contain any values.")
		return values

	if definition.mode in {"linspace", "range"}:
		parts = [part.strip() for part in spec.split(",") if part.strip()]
		if len(parts) != 3:
			raise SweepSpecError(
				f"Sweep '{definition.label}' expects 'start, stop, count' or 'start, stop, step'."
			)
		start = try_parse_float(parts[0])
		stop = try_parse_float(parts[1])
		third = try_parse_float(parts[2])
		if start is None or stop is None or third is None:
			raise SweepSpecError(f"Sweep '{definition.label}' requires numeric values.")

		if definition.mode == "linspace":
			count = int(round(third))
			if count <= 0:
				raise SweepSpecError(f"Sweep '{definition.label}' needs a positive sample count.")
			numbers = np.linspace(start, stop, count)
			return [format_numeric(float(value), template_value) for value in numbers]

		step = third
		if abs(step) < 1.0e-18:
			raise SweepSpecError(f"Sweep '{definition.label}' needs a non-zero step size.")
		values: list[str] = []
		current = start
		limit = 0
		if step > 0.0:
			while current <= stop + abs(step) * 1.0e-9:
				values.append(format_numeric(current, template_value))
				current += step
				limit += 1
				if limit > 100_000:
					raise SweepSpecError(f"Sweep '{definition.label}' expands to too many values.")
		else:
			while current >= stop - abs(step) * 1.0e-9:
				values.append(format_numeric(current, template_value))
				current += step
				limit += 1
				if limit > 100_000:
					raise SweepSpecError(f"Sweep '{definition.label}' expands to too many values.")
		return values

	raise SweepSpecError(f"Unsupported sweep mode: {definition.mode}")


def build_case_name(index: int, assignments: list[tuple[SweepDefinition, str]]) -> str:
	parts = [f"case_{index:04d}"]
	for definition, value in assignments:
		parts.append(f"{sanitize_token(definition.label, 24)}-{sanitize_token(value, 18)}")
	name = "__".join(parts)
	if len(name) > 180:
		name = name[:180].rstrip("-_.")
	return name


def expand_variants(definitions: list[SweepDefinition], entry_map: dict[str, CaseEntry]) -> list[VariantDescriptor]:
	enabled_definitions = [definition for definition in definitions if definition.enabled]
	if not enabled_definitions:
		return []

	sweep_values: list[list[str]] = []
	for definition in enabled_definitions:
		entry = entry_map[definition.entry_id]
		sweep_values.append(evaluate_sweep(definition, entry.value_text))

	variants: list[VariantDescriptor] = []
	for index, combination in enumerate(itertools.product(*sweep_values), start=1):
		assignments = list(zip(enabled_definitions, combination, strict=True))
		variants.append(
			VariantDescriptor(
				index=index,
				case_name=build_case_name(index, assignments),
				assignments=assignments,
			)
		)
	return variants


def apply_replacements(text: str, replacements: list[tuple[int, int, str]]) -> str:
	result = text
	for start, end, replacement in sorted(replacements, key=lambda item: item[0], reverse=True):
		result = result[:start] + replacement + result[end:]
	return result


class MatrixPreviewWidget(QWidget):
	def __init__(self, parent: QWidget | None = None) -> None:
		super().__init__(parent)
		self._pixmap: QPixmap | None = None

		layout = QVBoxLayout(self)
		layout.setContentsMargins(0, 0, 0, 0)

		self.info_label = QLabel("Select a numeric matrix asset to preview it.")
		self.info_label.setWordWrap(True)
		layout.addWidget(self.info_label)

		self.image_label = QLabel()
		self.image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
		self.image_label.setMinimumHeight(240)
		self.image_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
		layout.addWidget(self.image_label, stretch=1)

	def show_message(self, message: str) -> None:
		self.info_label.setText(message)
		self._pixmap = None
		self.image_label.clear()

	def show_file(self, path: Path) -> None:
		matrix = load_numeric_matrix(str(path))
		image = matrix_to_qimage(matrix)
		self._pixmap = QPixmap.fromImage(image)
		self.info_label.setText(f"{path.name}\n{describe_matrix(matrix)}")
		self._update_pixmap()

	def resizeEvent(self, event) -> None:  # type: ignore[override]
		super().resizeEvent(event)
		self._update_pixmap()

	def _update_pixmap(self) -> None:
		if self._pixmap is None:
			return
		scaled = self._pixmap.scaled(
			self.image_label.size(),
			Qt.AspectRatioMode.KeepAspectRatio,
			Qt.TransformationMode.SmoothTransformation,
		)
		self.image_label.setPixmap(scaled)


class SweepPreviewWidget(QOpenGLWidget):
	def __init__(self, parent: QWidget | None = None) -> None:
		super().__init__(parent)
		self.setMinimumHeight(220)
		self._points: list[tuple[float, float, float, tuple[float, float, float]]] = []
		self._axis_labels = ["x", "y", "z"]
		self._variant_count = 0
		self._azimuth = 35.0
		self._elevation = 25.0
		self._distance = 3.2
		self._last_mouse_pos: QPoint | None = None

	def set_variants(self, variants: list[VariantDescriptor]) -> None:
		self._variant_count = len(variants)
		if not variants:
			self._points = []
			self._axis_labels = ["x", "y", "z"]
			self.update()
			return

		labels = [assignment[0].label for assignment in variants[0].assignments]
		while len(labels) < 3:
			labels.append(f"axis-{len(labels) + 1}")
		self._axis_labels = labels[:3]

		dimensions: list[list[str]] = []
		for axis_index in range(min(3, len(labels))):
			values: list[str] = []
			for variant in variants:
				if axis_index < len(variant.assignments):
					values.append(variant.assignments[axis_index][1])
				else:
					values.append("0")
			dimensions.append(values)

		axis_coordinates = [self._normalize_dimension(values) for values in dimensions]
		self._points = []
		max_count = max(1, len(variants) - 1)
		for index, variant in enumerate(variants):
			x_value = axis_coordinates[0][index] if axis_coordinates else 0.5
			y_value = axis_coordinates[1][index] if len(axis_coordinates) > 1 else 0.5
			z_value = axis_coordinates[2][index] if len(axis_coordinates) > 2 else 0.5
			hue = index / max_count
			color = self._color_from_hue(hue)
			self._points.append((x_value, y_value, z_value, color))
		self.update()

	def _normalize_dimension(self, values: list[str]) -> list[float]:
		numeric = [try_parse_float(value) for value in values]
		if all(value is not None for value in numeric):
			floats = [float(value) for value in numeric if value is not None]
			min_value = min(floats)
			max_value = max(floats)
			if abs(max_value - min_value) < 1.0e-12:
				return [0.5 for _ in floats]
			return [(value - min_value) / (max_value - min_value) for value in floats]

		ordered = sorted(dict.fromkeys(values))
		denominator = max(1, len(ordered) - 1)
		mapping = {value: index / denominator for index, value in enumerate(ordered)}
		return [mapping[value] for value in values]

	def _color_from_hue(self, hue: float) -> tuple[float, float, float]:
		angle = hue * 2.0 * math.pi
		return (
			0.5 + 0.5 * math.sin(angle),
			0.5 + 0.5 * math.sin(angle + 2.0 * math.pi / 3.0),
			0.5 + 0.5 * math.sin(angle + 4.0 * math.pi / 3.0),
		)

	def initializeGL(self) -> None:  # type: ignore[override]
		if not HAVE_PYOPENGL:
			return
		GL.glEnable(GL.GL_DEPTH_TEST)
		GL.glEnable(GL.GL_POINT_SMOOTH)
		GL.glPointSize(7.0)

	def resizeGL(self, width: int, height: int) -> None:  # type: ignore[override]
		if HAVE_PYOPENGL:
			GL.glViewport(0, 0, width, max(1, height))

	def paintGL(self) -> None:  # type: ignore[override]
		if not HAVE_PYOPENGL:
			painter = QPainter(self)
			painter.fillRect(self.rect(), QColor("#111827"))
			painter.setPen(QColor("#e5e7eb"))
			painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "PyOpenGL not available")
			painter.end()
			return

		GL.glClearColor(0.07, 0.09, 0.14, 1.0)
		GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

		aspect = max(1.0e-6, self.width() / max(1.0, float(self.height())))
		GL.glMatrixMode(GL.GL_PROJECTION)
		GL.glLoadIdentity()
		GLU.gluPerspective(45.0, aspect, 0.1, 100.0)

		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glLoadIdentity()
		GL.glTranslatef(0.0, 0.0, -self._distance)
		GL.glRotatef(self._elevation, 1.0, 0.0, 0.0)
		GL.glRotatef(self._azimuth, 0.0, 1.0, 0.0)
		GL.glTranslatef(-0.5, -0.5, -0.5)

		self._draw_axes()
		self._draw_points()

		painter = QPainter(self)
		painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
		painter.setPen(QColor("#f9fafb"))
		painter.drawText(12, 20, f"Variants: {self._variant_count}")
		painter.drawText(12, 40, f"Axes: {self._axis_labels[0]}, {self._axis_labels[1]}, {self._axis_labels[2]}")
		painter.drawText(12, self.height() - 12, "Drag to orbit, mouse wheel to zoom")
		painter.end()

	def _draw_axes(self) -> None:
		GL.glLineWidth(2.0)
		GL.glBegin(GL.GL_LINES)
		GL.glColor3f(0.9, 0.3, 0.3)
		GL.glVertex3f(0.0, 0.0, 0.0)
		GL.glVertex3f(1.0, 0.0, 0.0)
		GL.glColor3f(0.3, 0.9, 0.3)
		GL.glVertex3f(0.0, 0.0, 0.0)
		GL.glVertex3f(0.0, 1.0, 0.0)
		GL.glColor3f(0.3, 0.6, 1.0)
		GL.glVertex3f(0.0, 0.0, 0.0)
		GL.glVertex3f(0.0, 0.0, 1.0)

		GL.glColor4f(0.4, 0.4, 0.4, 1.0)
		for corner in (0.0, 1.0):
			GL.glVertex3f(corner, 0.0, 0.0)
			GL.glVertex3f(corner, 1.0, 0.0)
			GL.glVertex3f(0.0, corner, 0.0)
			GL.glVertex3f(1.0, corner, 0.0)
			GL.glVertex3f(0.0, 0.0, corner)
			GL.glVertex3f(1.0, 0.0, corner)
		GL.glEnd()

	def _draw_points(self) -> None:
		if not self._points:
			return
		GL.glBegin(GL.GL_POINTS)
		for x_value, y_value, z_value, color in self._points:
			GL.glColor3f(*color)
			GL.glVertex3f(x_value, y_value, z_value)
		GL.glEnd()

	def mousePressEvent(self, event) -> None:  # type: ignore[override]
		if event.button() == Qt.MouseButton.LeftButton:
			self._last_mouse_pos = event.position().toPoint()

	def mouseMoveEvent(self, event) -> None:  # type: ignore[override]
		if self._last_mouse_pos is None:
			return
		position = event.position().toPoint()
		delta = position - self._last_mouse_pos
		self._azimuth += delta.x() * 0.6
		self._elevation += delta.y() * 0.6
		self._elevation = max(-89.0, min(89.0, self._elevation))
		self._last_mouse_pos = position
		self.update()

	def mouseReleaseEvent(self, event) -> None:  # type: ignore[override]
		self._last_mouse_pos = None
		super().mouseReleaseEvent(event)

	def wheelEvent(self, event) -> None:  # type: ignore[override]
		delta = event.angleDelta().y() / 120.0
		self._distance = max(1.4, min(10.0, self._distance - delta * 0.2))
		self.update()


class SpawnMainWindow(QMainWindow):
	def __init__(self, case_root: Path | None = None) -> None:
		super().__init__()
		self.case_root = case_root or default_case_root()
		self.output_root = self.case_root.parent / f"{self.case_root.name}_spawned"
		self.parsed_files: dict[str, ParsedFile] = {}
		self.entry_by_id: dict[str, CaseEntry] = {}
		self.current_file_key: str | None = None
		self._sweep_table_updating = False

		self.setWindowTitle("laserbeamFoam Spawn")
		self.resize(1680, 980)

		self._build_ui()
		self._build_toolbar()
		self.setStatusBar(QStatusBar())

		if is_case_directory(self.case_root):
			self.load_case(self.case_root)
		else:
			self.statusBar().showMessage("Choose a laserbeamFoam case directory to begin.")

	def _build_toolbar(self) -> None:
		toolbar = QToolBar("Main")
		toolbar.setMovable(False)
		self.addToolBar(toolbar)

		open_action = QAction("Open Case", self)
		open_action.triggered.connect(self.browse_case)
		toolbar.addAction(open_action)

		refresh_action = QAction("Refresh", self)
		refresh_action.triggered.connect(self.reload_case)
		toolbar.addAction(refresh_action)

		generate_action = QAction("Generate", self)
		generate_action.triggered.connect(self.generate_cases)
		toolbar.addAction(generate_action)

	def _build_ui(self) -> None:
		container = QWidget()
		self.setCentralWidget(container)
		outer_layout = QVBoxLayout(container)
		outer_layout.setContentsMargins(8, 8, 8, 8)

		root_form = QFormLayout()
		self.case_path_edit = QLineEdit(str(self.case_root))
		self.case_path_edit.setReadOnly(True)
		case_button = QPushButton("Browse...")
		case_button.clicked.connect(self.browse_case)
		case_row = QHBoxLayout()
		case_row.addWidget(self.case_path_edit, stretch=1)
		case_row.addWidget(case_button)
		root_form.addRow("Template Case", wrap_layout(case_row))

		self.output_path_edit = QLineEdit(str(self.output_root))
		output_button = QPushButton("Browse...")
		output_button.clicked.connect(self.browse_output_root)
		output_row = QHBoxLayout()
		output_row.addWidget(self.output_path_edit, stretch=1)
		output_row.addWidget(output_button)
		root_form.addRow("Output Root", wrap_layout(output_row))
		outer_layout.addLayout(root_form)

		root_splitter = QSplitter(Qt.Orientation.Horizontal)
		outer_layout.addWidget(root_splitter, stretch=1)

		left_panel = QWidget()
		left_layout = QVBoxLayout(left_panel)
		left_layout.setContentsMargins(0, 0, 0, 0)

		file_group = QGroupBox("Case Files")
		file_layout = QVBoxLayout(file_group)
		self.file_tree = QTreeWidget()
		self.file_tree.setHeaderLabels(["File", "Kind"])
		self.file_tree.setAlternatingRowColors(True)
		self.file_tree.itemSelectionChanged.connect(self.on_file_selection_changed)
		file_layout.addWidget(self.file_tree)
		left_layout.addWidget(file_group, stretch=1)

		root_splitter.addWidget(left_panel)

		right_splitter = QSplitter(Qt.Orientation.Vertical)
		root_splitter.addWidget(right_splitter)
		root_splitter.setStretchFactor(1, 1)

		entries_group = QGroupBox("Editable Entries")
		entries_layout = QVBoxLayout(entries_group)
		entry_toolbar = QHBoxLayout()
		self.entry_summary_label = QLabel("Select a parsed OpenFOAM input file.")
		entry_toolbar.addWidget(self.entry_summary_label, stretch=1)
		add_button = QPushButton("Add Selected To Sweeps")
		add_button.clicked.connect(self.add_selected_entries_to_sweeps)
		entry_toolbar.addWidget(add_button)
		entries_layout.addLayout(entry_toolbar)

		self.entry_table = QTableWidget(0, 2)
		self.entry_table.setHorizontalHeaderLabels(["Entry", "Value"])
		self.entry_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
		self.entry_table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
		self.entry_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
		self.entry_table.verticalHeader().setVisible(False)
		self.entry_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
		self.entry_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
		entries_layout.addWidget(self.entry_table, stretch=1)
		right_splitter.addWidget(entries_group)

		preview_splitter = QSplitter(Qt.Orientation.Horizontal)
		right_splitter.addWidget(preview_splitter)

		raw_group = QGroupBox("Raw Preview")
		raw_layout = QVBoxLayout(raw_group)
		self.raw_preview = QPlainTextEdit()
		self.raw_preview.setReadOnly(True)
		raw_layout.addWidget(self.raw_preview)
		preview_splitter.addWidget(raw_group)

		side_preview = QWidget()
		side_preview_layout = QVBoxLayout(side_preview)
		side_preview_layout.setContentsMargins(0, 0, 0, 0)

		numeric_group = QGroupBox("Numeric Asset Preview")
		numeric_layout = QVBoxLayout(numeric_group)
		self.matrix_preview = MatrixPreviewWidget()
		numeric_layout.addWidget(self.matrix_preview)
		side_preview_layout.addWidget(numeric_group, stretch=1)

		log_group = QGroupBox("Log Summary")
		log_layout = QVBoxLayout(log_group)
		self.log_preview = QPlainTextEdit()
		self.log_preview.setReadOnly(True)
		log_layout.addWidget(self.log_preview)
		side_preview_layout.addWidget(log_group, stretch=1)
		preview_splitter.addWidget(side_preview)
		preview_splitter.setStretchFactor(0, 1)

		bottom_splitter = QSplitter(Qt.Orientation.Horizontal)
		right_splitter.addWidget(bottom_splitter)

		sweep_group = QGroupBox("Sweep Definitions")
		sweep_layout = QVBoxLayout(sweep_group)
		sweep_toolbar = QHBoxLayout()
		remove_button = QPushButton("Remove Selected")
		remove_button.clicked.connect(self.remove_selected_sweeps)
		clear_button = QPushButton("Clear")
		clear_button.clicked.connect(self.clear_sweeps)
		sweep_toolbar.addWidget(remove_button)
		sweep_toolbar.addWidget(clear_button)
		sweep_toolbar.addStretch(1)
		sweep_layout.addLayout(sweep_toolbar)

		self.sweep_table = QTableWidget(0, 6)
		self.sweep_table.setHorizontalHeaderLabels(["On", "Entry", "Current", "Mode", "Spec", "Label"])
		self.sweep_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
		self.sweep_table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
		self.sweep_table.verticalHeader().setVisible(False)
		self.sweep_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
		self.sweep_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
		self.sweep_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
		self.sweep_table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
		self.sweep_table.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeMode.Stretch)
		self.sweep_table.horizontalHeader().setSectionResizeMode(5, QHeaderView.ResizeMode.ResizeToContents)
		self.sweep_table.itemChanged.connect(self.update_variants)
		sweep_layout.addWidget(self.sweep_table, stretch=1)

		self.variant_summary_label = QLabel("No sweeps defined.")
		sweep_layout.addWidget(self.variant_summary_label)

		self.variant_table = QTableWidget(0, 0)
		self.variant_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
		self.variant_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
		self.variant_table.verticalHeader().setVisible(False)
		sweep_layout.addWidget(self.variant_table, stretch=1)
		bottom_splitter.addWidget(sweep_group)

		preview_group = QGroupBox("Array Preview")
		preview_layout = QVBoxLayout(preview_group)
		self.opengl_preview = SweepPreviewWidget()
		preview_layout.addWidget(self.opengl_preview)
		bottom_splitter.addWidget(preview_group)
		bottom_splitter.setStretchFactor(0, 1)

		right_splitter.setStretchFactor(0, 3)
		right_splitter.setStretchFactor(1, 2)
		right_splitter.setStretchFactor(2, 3)

	def browse_case(self) -> None:
		chosen = QFileDialog.getExistingDirectory(self, "Choose Template Case", str(self.case_root))
		if not chosen:
			return
		self.load_case(Path(chosen).resolve())

	def browse_output_root(self) -> None:
		chosen = QFileDialog.getExistingDirectory(self, "Choose Output Root", self.output_path_edit.text())
		if not chosen:
			return
		self.output_root = Path(chosen).resolve()
		self.output_path_edit.setText(str(self.output_root))

	def reload_case(self) -> None:
		self.load_case(self.case_root)

	def load_case(self, case_root: Path) -> None:
		if not is_case_directory(case_root):
			QMessageBox.warning(self, "Invalid Case", f"{case_root} does not look like a laserbeamFoam case.")
			return

		QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
		try:
			self.case_root = case_root.resolve()
			self.output_root = self.case_root.parent / f"{self.case_root.name}_spawned"
			self.case_path_edit.setText(str(self.case_root))
			self.output_path_edit.setText(str(self.output_root))
			self.parsed_files = scan_case(self.case_root)
			self.entry_by_id = {
				entry.identifier: entry
				for parsed_file in self.parsed_files.values()
				for entry in parsed_file.entries
			}
			self.populate_file_tree()
			self.clear_previews()
			self.clear_sweeps()
			self.statusBar().showMessage(
				f"Loaded {len(self.parsed_files)} files and {len(self.entry_by_id)} editable entries from {self.case_root.name}."
			)
		finally:
			QApplication.restoreOverrideCursor()

	def populate_file_tree(self) -> None:
		self.file_tree.clear()
		nodes: dict[tuple[str, ...], QTreeWidgetItem] = {(): self.file_tree.invisibleRootItem()}

		for key in sorted(self.parsed_files):
			parsed_file = self.parsed_files[key]
			relative_parts = parsed_file.relative_path.parts
			parent_key: tuple[str, ...] = ()
			for depth, part in enumerate(relative_parts):
				node_key = relative_parts[: depth + 1]
				if node_key in nodes:
					parent_key = node_key
					continue

				if depth == len(relative_parts) - 1:
					item = QTreeWidgetItem([part, parsed_file.kind])
					item.setData(0, Qt.ItemDataRole.UserRole, key)
				else:
					item = QTreeWidgetItem([part, "dir"])
				nodes[parent_key].addChild(item)
				nodes[node_key] = item
				parent_key = node_key

		self.file_tree.expandToDepth(1)

		for index in range(self.file_tree.topLevelItemCount()):
			top = self.file_tree.topLevelItem(index)
			if top.childCount() > 0:
				candidate = self._first_file_item(top)
				if candidate is not None:
					self.file_tree.setCurrentItem(candidate)
					break

	def _first_file_item(self, item: QTreeWidgetItem) -> QTreeWidgetItem | None:
		key = item.data(0, Qt.ItemDataRole.UserRole)
		if isinstance(key, str):
			return item
		for child_index in range(item.childCount()):
			child = item.child(child_index)
			candidate = self._first_file_item(child)
			if candidate is not None:
				return candidate
		return None

	def clear_previews(self) -> None:
		self.current_file_key = None
		self.entry_summary_label.setText("Select a parsed OpenFOAM input file.")
		self.entry_table.setRowCount(0)
		self.raw_preview.clear()
		self.matrix_preview.show_message("Select a numeric matrix asset to preview it.")
		self.log_preview.clear()

	def on_file_selection_changed(self) -> None:
		selected = self.file_tree.selectedItems()
		if not selected:
			return
		key = selected[0].data(0, Qt.ItemDataRole.UserRole)
		if not isinstance(key, str):
			return
		self.show_file(key)

	def show_file(self, key: str) -> None:
		parsed_file = self.parsed_files[key]
		self.current_file_key = key

		entry_count = len(parsed_file.entries)
		parse_note = f" | parser warning: {parsed_file.parse_error}" if parsed_file.parse_error else ""
		self.entry_summary_label.setText(
			f"{parsed_file.relative_path.as_posix()} | kind={parsed_file.kind} | entries={entry_count}{parse_note}"
		)

		self.populate_entries(parsed_file.entries)

		if parsed_file.text is not None:
			self.raw_preview.setPlainText(parsed_file.text)
		elif parsed_file.kind == "log":
			self.raw_preview.setPlainText(parsed_file.path.read_text(encoding="utf-8", errors="ignore"))
		else:
			self.raw_preview.setPlainText(parsed_file.path.read_text(encoding="utf-8", errors="ignore"))

		if parsed_file.kind == "numeric":
			self.matrix_preview.show_file(parsed_file.path)
		else:
			self.matrix_preview.show_message("Selected file is not a numeric matrix asset.")

		if parsed_file.kind == "log":
			try:
				self.log_preview.setPlainText(summarize_log(parsed_file.path))
			except Exception as exc:
				self.log_preview.setPlainText(f"Could not summarize {parsed_file.path.name}:\n{exc}")
		else:
			self.log_preview.setPlainText("Select a log.* file to reuse the shared plot parser summary.")

	def populate_entries(self, entries: list[CaseEntry]) -> None:
		self.entry_table.setRowCount(len(entries))
		for row, entry in enumerate(entries):
			path_item = QTableWidgetItem(entry.path_label)
			path_item.setData(Qt.ItemDataRole.UserRole, entry.identifier)
			path_item.setToolTip(entry.identifier)
			value_item = QTableWidgetItem(entry.display_value)
			value_item.setToolTip(entry.value_text)
			self.entry_table.setItem(row, 0, path_item)
			self.entry_table.setItem(row, 1, value_item)

	def add_selected_entries_to_sweeps(self) -> None:
		selected_rows = sorted({index.row() for index in self.entry_table.selectionModel().selectedRows()})
		if not selected_rows:
			QMessageBox.information(self, "No Selection", "Select one or more entries first.")
			return

		existing_ids = {self.sweep_table.item(row, 1).data(Qt.ItemDataRole.UserRole) for row in range(self.sweep_table.rowCount())}
		for row in selected_rows:
			entry_id = self.entry_table.item(row, 0).data(Qt.ItemDataRole.UserRole)
			if entry_id in existing_ids:
				continue
			entry = self.entry_by_id[entry_id]
			self._append_sweep_row(entry)

		self.update_variants()

	def _append_sweep_row(self, entry: CaseEntry) -> None:
		self._sweep_table_updating = True
		try:
			row = self.sweep_table.rowCount()
			self.sweep_table.insertRow(row)

			enabled_item = QTableWidgetItem()
			enabled_item.setFlags(enabled_item.flags() | Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
			enabled_item.setCheckState(Qt.CheckState.Checked)
			self.sweep_table.setItem(row, 0, enabled_item)

			entry_item = QTableWidgetItem(entry.identifier)
			entry_item.setData(Qt.ItemDataRole.UserRole, entry.identifier)
			entry_item.setToolTip(entry.identifier)
			self.sweep_table.setItem(row, 1, entry_item)

			current_item = QTableWidgetItem(entry.display_value)
			current_item.setToolTip(entry.value_text)
			current_item.setFlags(current_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
			self.sweep_table.setItem(row, 2, current_item)

			mode_combo = QComboBox()
			mode_combo.addItems(["list", "linspace", "range"])
			mode_combo.currentTextChanged.connect(self.update_variants)
			self.sweep_table.setCellWidget(row, 3, mode_combo)

			spec_item = QTableWidgetItem(entry.value_text.strip())
			self.sweep_table.setItem(row, 4, spec_item)

			label_item = QTableWidgetItem(sanitize_token(entry.short_label, 24))
			self.sweep_table.setItem(row, 5, label_item)
		finally:
			self._sweep_table_updating = False

	def remove_selected_sweeps(self) -> None:
		rows = sorted({index.row() for index in self.sweep_table.selectionModel().selectedRows()}, reverse=True)
		for row in rows:
			self.sweep_table.removeRow(row)
		self.update_variants()

	def clear_sweeps(self) -> None:
		self._sweep_table_updating = True
		try:
			self.sweep_table.setRowCount(0)
		finally:
			self._sweep_table_updating = False
		self.variant_table.setRowCount(0)
		self.variant_table.setColumnCount(0)
		self.variant_summary_label.setText("No sweeps defined.")
		self.opengl_preview.set_variants([])

	def collect_sweeps(self) -> list[SweepDefinition]:
		definitions: list[SweepDefinition] = []
		for row in range(self.sweep_table.rowCount()):
			enabled_item = self.sweep_table.item(row, 0)
			entry_item = self.sweep_table.item(row, 1)
			spec_item = self.sweep_table.item(row, 4)
			label_item = self.sweep_table.item(row, 5)
			mode_combo = self.sweep_table.cellWidget(row, 3)

			if not isinstance(mode_combo, QComboBox):
				continue
			if enabled_item is None or entry_item is None or spec_item is None or label_item is None:
				continue

			definitions.append(
				SweepDefinition(
					entry_id=entry_item.data(Qt.ItemDataRole.UserRole),
					label=label_item.text().strip() or sanitize_token(entry_item.text().split("::", maxsplit=1)[-1], 24),
					mode=mode_combo.currentText(),
					spec=spec_item.text(),
					enabled=enabled_item.checkState() == Qt.CheckState.Checked,
				)
			)
		return definitions

	def update_variants(self) -> None:
		if self._sweep_table_updating:
			return

		definitions = self.collect_sweeps()
		if not definitions:
			self.variant_table.setRowCount(0)
			self.variant_table.setColumnCount(0)
			self.variant_summary_label.setText("No sweeps defined.")
			self.opengl_preview.set_variants([])
			return

		try:
			variants = expand_variants(definitions, self.entry_by_id)
		except SweepSpecError as exc:
			self.variant_summary_label.setText(str(exc))
			self.variant_table.setRowCount(0)
			self.variant_table.setColumnCount(0)
			self.opengl_preview.set_variants([])
			return
		except KeyError as exc:
			self.variant_summary_label.setText(f"Unknown entry in sweep table: {exc}")
			self.variant_table.setRowCount(0)
			self.variant_table.setColumnCount(0)
			self.opengl_preview.set_variants([])
			return

		self.variant_summary_label.setText(f"Previewing {len(variants)} variants.")
		self.populate_variant_table(variants)
		self.opengl_preview.set_variants(variants)

	def populate_variant_table(self, variants: list[VariantDescriptor]) -> None:
		if not variants:
			self.variant_table.setRowCount(0)
			self.variant_table.setColumnCount(0)
			return

		headers = ["Case"] + [assignment[0].label for assignment in variants[0].assignments]
		preview_rows = min(len(variants), 256)
		self.variant_table.setColumnCount(len(headers))
		self.variant_table.setHorizontalHeaderLabels(headers)
		self.variant_table.setRowCount(preview_rows)

		for row in range(preview_rows):
			variant = variants[row]
			self.variant_table.setItem(row, 0, QTableWidgetItem(variant.case_name))
			for column, (_definition, value) in enumerate(variant.assignments, start=1):
				self.variant_table.setItem(row, column, QTableWidgetItem(value))

		self.variant_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
		for column in range(1, len(headers)):
			self.variant_table.horizontalHeader().setSectionResizeMode(column, QHeaderView.ResizeMode.ResizeToContents)

		if len(variants) > preview_rows:
			self.variant_summary_label.setText(
				f"Previewing {preview_rows} of {len(variants)} variants. Generation still uses the full set."
			)

	def generate_cases(self) -> None:
		definitions = self.collect_sweeps()
		if not definitions:
			QMessageBox.information(self, "No Sweeps", "Add at least one sweep entry before generating cases.")
			return

		try:
			variants = expand_variants(definitions, self.entry_by_id)
		except SweepSpecError as exc:
			QMessageBox.warning(self, "Invalid Sweep", str(exc))
			return

		if not variants:
			QMessageBox.information(self, "No Variants", "No enabled sweep definitions were found.")
			return

		output_root = Path(self.output_path_edit.text()).expanduser().resolve()
		output_root.mkdir(parents=True, exist_ok=True)

		progress = QProgressDialog("Generating cases...", "Abort", 0, len(variants), self)
		progress.setWindowTitle("Generating")
		progress.setWindowModality(Qt.WindowModality.WindowModal)
		progress.show()

		manifest = {
			"template_case": str(self.case_root),
			"output_root": str(output_root),
			"generated_at": datetime.now().isoformat(timespec="seconds"),
			"sweeps": [definition.__dict__ for definition in definitions],
			"variants": [],
		}

		try:
			for index, variant in enumerate(variants, start=1):
				if progress.wasCanceled():
					break

				case_dir = output_root / variant.case_name
				if case_dir.exists():
					raise FileExistsError(f"Refusing to overwrite existing directory: {case_dir}")

				shutil.copytree(self.case_root, case_dir, ignore=self._build_ignore_function())

				replacements_by_file: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
				variant_manifest = {
					"case_name": variant.case_name,
					"index": variant.index,
					"assignments": {},
				}

				for definition, value in variant.assignments:
					entry = self.entry_by_id[definition.entry_id]
					replacements_by_file[entry.relative_path.as_posix()].append((entry.value_start, entry.value_end, value))
					variant_manifest["assignments"][definition.entry_id] = value

				for relative_text, replacements in replacements_by_file.items():
					target = case_dir / relative_text
					original_text = target.read_text(encoding="utf-8", errors="ignore")
					updated_text = apply_replacements(original_text, replacements)
					target.write_text(updated_text, encoding="utf-8")

				(case_dir / "spawn_variant.json").write_text(
					json.dumps(variant_manifest, indent=2, sort_keys=True),
					encoding="utf-8",
				)
				manifest["variants"].append(variant_manifest)
				progress.setValue(index)
				progress.setLabelText(f"Generated {variant.case_name}")

			(output_root / "spawn_manifest.json").write_text(
				json.dumps(manifest, indent=2, sort_keys=True),
				encoding="utf-8",
			)
		except Exception as exc:
			QMessageBox.critical(self, "Generation Failed", str(exc))
			return
		finally:
			progress.close()

		self.statusBar().showMessage(f"Generated {len(manifest['variants'])} cases under {output_root}")
		QMessageBox.information(
			self,
			"Generation Complete",
			f"Generated {len(manifest['variants'])} cases under\n{output_root}\n\nA manifest was written to spawn_manifest.json.",
		)

	def _build_ignore_function(self):
		def _ignore(directory: str, names: list[str]) -> list[str]:
			ignored: list[str] = []
			dir_path = Path(directory)
			for name in names:
				full_path = dir_path / name
				try:
					relative_path = full_path.relative_to(self.case_root)
				except ValueError:
					continue
				if should_skip_generated_path(relative_path):
					ignored.append(name)
			return ignored

		return _ignore


def wrap_layout(layout: QHBoxLayout) -> QWidget:
	widget = QWidget()
	widget.setLayout(layout)
	return widget


def configure_surface_format() -> None:
	fmt = QSurfaceFormat()
	fmt.setDepthBufferSize(24)
	fmt.setSamples(4)
	QSurfaceFormat.setDefaultFormat(fmt)


def main(argv: list[str]) -> int:
	configure_surface_format()
	app = QApplication(argv)
	app.setApplicationName("laserbeamFoam Spawn")
	app.setOrganizationName("laserbeamFoam")

	case_root: Path | None = None
	if len(argv) > 1:
		case_root = Path(argv[1]).expanduser().resolve()

	window = SpawnMainWindow(case_root)
	window.show()
	return app.exec()


if __name__ == "__main__":
	raise SystemExit(main(sys.argv))
