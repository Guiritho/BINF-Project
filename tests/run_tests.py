import argparse
import contextlib
import csv
import enum
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_SRC_DIR = PROJECT_ROOT / "src"
DEFAULT_TESTS_DIR = PROJECT_ROOT / "data" / "tests" / "230626_current"

SUBPROCESS_TIMEOUT_SECONDS = 300
CLEANED_TOLERANCE = 1e-4
TOM_TOLERANCE = 1e-3
EIGEN_TOLERANCE = 1e-2
PROGRAM_NOT_FOUND_MESSAGE = "program not found in src/"

COLOR_RESET = "\033[0m"
COLOR_BOLD = "\033[1m"
COLOR_GREEN = "\033[32m"
COLOR_RED = "\033[31m"
COLOR_YELLOW = "\033[33m"
COLOR_CYAN = "\033[36m"
COLOR_GRAY = "\033[90m"
CLEAR_LINE = "\033[K"

STATUS_COLORS = {
    "PASS": COLOR_GREEN,
    "FAIL": COLOR_RED,
    "ERROR": COLOR_RED,
    "SKIP": COLOR_YELLOW,
}

STATUS_LABEL_WIDTH = 5
PROGRESS_BAR_WIDTH = 30
FALLBACK_TERMINAL_WIDTH = 100
TRUNCATION_SUFFIX = "..."

DATASETS = ["small", "medium", "large", "full"]
GENE_COUNTS = [50, 100, 250, 500, 750, 1000]
STATIC_L_VALUES = [0.6, 0.9]
TAU_VALUES = [1, 15, 30]
ADAPTIVE_AND_DYNAMIC_OPTIONS = ["adaptive", "dynamic"]

EXPECTED_BETA = {
    ("small", 50): 13, ("small", 100): 11, ("small", 250): 12,
    ("small", 500): 10, ("small", 750): 11, ("small", 1000): 10,
    ("medium", 50): 1, ("medium", 100): 10, ("medium", 250): 10,
    ("medium", 500): 9, ("medium", 750): 9, ("medium", 1000): 10,
    ("large", 50): 10, ("large", 100): 10, ("large", 250): 9,
    ("large", 500): 9, ("large", 750): 9, ("large", 1000): 8,
    ("full", 50): 2, ("full", 100): 9, ("full", 250): 8,
    ("full", 500): 9, ("full", 750): 8, ("full", 1000): 8,
}


class TestStatus(enum.Enum):
    PASSED = "PASS"
    FAILED = "FAIL"
    ERROR = "ERROR"
    SKIPPED = "SKIP"


@dataclass
class TestOutcome:
    name: str
    status: TestStatus
    message: str
    duration: float = 0.0


@dataclass
class TestContext:
    src_dir: Path
    tests_dir: Path
    timeout: int
    verbose: bool
    use_color: bool
    is_tty: bool


# --- Generic helpers --------------------------------------------------------

@contextlib.contextmanager
def staged_directory(tests_dir, input_filenames):
    """Create a temp directory containing copies of the given reference files."""
    with tempfile.TemporaryDirectory(prefix="wgcna_test_") as tmpdir:
        workdir = Path(tmpdir)
        for filename in input_filenames:
            shutil.copy(tests_dir / filename, workdir / filename)
        yield workdir


def run_program(context, program_name, args):
    script_path = context.src_dir / program_name
    if not script_path.exists():
        return None, PROGRAM_NOT_FOUND_MESSAGE

    command = ["python3", str(script_path)] + args
    result = subprocess.run(
        command, capture_output=True, text=True, timeout=context.timeout,
    )
    return result, None


def run_and_check(context, program_name, args, test_name):
    """Run a program, returning either an error outcome or the completed result."""
    result, error = run_program(context, program_name, args)
    if error == PROGRAM_NOT_FOUND_MESSAGE:
        return None, TestOutcome(test_name, TestStatus.SKIPPED, error)
    if error is not None:
        return None, TestOutcome(test_name, TestStatus.ERROR, error)
    if result.returncode != 0:
        return None, TestOutcome(test_name, TestStatus.ERROR, result.stderr.strip())

    return result, None


def timed(check_fn, *args):
    """Call check_fn(*args), returning its TestOutcome with duration filled in."""
    start = time.perf_counter()
    outcome = check_fn(*args)
    outcome.duration = time.perf_counter() - start

    return outcome


# --- File readers ------------------------------------------------------------

def read_numeric_matrix(filepath):
    """Read a CSV with no header and no row ids (TOM files)."""
    return np.loadtxt(filepath, delimiter=",")


def read_cleaned_csv(filepath):
    """Read a cleaned-dataset CSV (header + row ids), return data values only."""
    with open(filepath) as handle:
        reader = csv.reader(handle)
        next(reader)
        rows = [row[1:] for row in reader]

    return np.array(rows, dtype=float)


def read_groups_csv(filepath):
    with open(filepath) as handle:
        return [int(line.strip()) for line in handle if line.strip()]


def read_eigen_csv(filepath):
    """Return a list of (group_id, values) pairs, one per line."""
    rows = []
    with open(filepath) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            group_id = int(parts[0])
            values = np.array(parts[1:], dtype=float)
            rows.append((group_id, values))

    return rows


# --- Comparators --------------------------------------------------------------

def compare_matrices(produced, reference, tolerance):
    if produced.shape != reference.shape:
        return False, f"shape mismatch: {produced.shape} vs {reference.shape}"

    max_diff = float(np.max(np.abs(produced - reference)))
    is_equal = max_diff <= tolerance

    return is_equal, f"max diff {max_diff:.6g} (tolerance {tolerance:g})"


def is_same_partition(produced, reference):
    """Check whether two group label sequences describe the same partition."""
    mapping = {}
    for produced_id, reference_id in zip(produced, reference):
        if produced_id in mapping:
            if mapping[produced_id] != reference_id:
                return False
        else:
            mapping[produced_id] = reference_id

    return len(set(mapping.values())) == len(set(mapping.keys()))


def compare_groups(produced, reference):
    if len(produced) != len(reference):
        return False, f"length mismatch: {len(produced)} vs {len(reference)}"
    if produced == reference:
        return True, "exact match"
    if is_same_partition(produced, reference):
        return True, "same partition up to relabeling"

    return False, "different partitions"


def compare_eigen(produced, reference, tolerance):
    if len(produced) != len(reference):
        return False, f"group count mismatch: {len(produced)} vs {len(reference)}"

    max_diff = 0.0
    for (produced_id, produced_values), (reference_id, reference_values) in zip(produced, reference):
        if produced_id != reference_id:
            return False, f"group id mismatch: {produced_id} vs {reference_id}"
        if produced_values.shape != reference_values.shape:
            return False, f"group {produced_id}: length mismatch"

        diff_same_sign = np.max(np.abs(produced_values - reference_values))
        diff_flipped_sign = np.max(np.abs(produced_values + reference_values))
        max_diff = max(max_diff, min(diff_same_sign, diff_flipped_sign))

    is_equal = max_diff <= tolerance
    return is_equal, f"max diff {max_diff:.6g} (tolerance {tolerance:g}, sign-agnostic)"


# --- Palier 1: clean_dataset --------------------------------------------------

def check_clean_dataset(context, dataset, gene_count):
    dataname = f"mrna2_{dataset}.csv"
    test_name = f"clean_dataset {dataname} K={gene_count}"

    with staged_directory(context.tests_dir, [dataname]) as workdir:
        args = [str(workdir), dataname, str(gene_count)]
        result, outcome = run_and_check(context, "clean_dataset.py", args, test_name)
        if outcome is not None:
            return outcome

        output_filename = f"{dataname}_cleaned_{gene_count}.csv"
        produced_path = workdir / output_filename
        if not produced_path.exists():
            return TestOutcome(test_name, TestStatus.FAILED, f"missing output {output_filename}")

        produced = read_cleaned_csv(produced_path)
        reference = read_cleaned_csv(context.tests_dir / output_filename)
        is_equal, message = compare_matrices(produced, reference, CLEANED_TOLERANCE)

        status = TestStatus.PASSED if is_equal else TestStatus.FAILED
        return TestOutcome(test_name, status, message)


def run_clean_dataset_tests(context, dataset_params):
    for dataset, gene_count, _ in dataset_params:
        yield timed(check_clean_dataset, context, dataset, gene_count)


def count_clean_dataset_tests(dataset_params):
    return len(dataset_params)


# --- Palier 2: build_graph (POW + TOM) ----------------------------------------

def check_build_graph_pow(context, dataset, gene_count, expected_beta):
    dataname = f"mrna2_{dataset}.csv_cleaned_{gene_count}.csv"
    test_name = f"build_graph POW {dataname}"

    with staged_directory(context.tests_dir, [dataname]) as workdir:
        args = [str(workdir), dataname, "TOM"]
        result, outcome = run_and_check(context, "build_graph.py", args, test_name)
        if outcome is not None:
            return outcome

        produced_beta = result.stdout.strip()
        if produced_beta != str(expected_beta):
            return TestOutcome(test_name, TestStatus.FAILED, f"beta={produced_beta}, expected {expected_beta}")

        return TestOutcome(test_name, TestStatus.PASSED, f"beta={expected_beta}")


def check_build_graph_tom(context, dataset, gene_count, beta):
    dataname = f"mrna2_{dataset}.csv_cleaned_{gene_count}.csv"
    test_name = f"build_graph TOM {dataname} beta={beta}"

    with staged_directory(context.tests_dir, [dataname]) as workdir:
        args = [str(workdir), dataname, "TOM", str(beta)]
        result, outcome = run_and_check(context, "build_graph.py", args, test_name)
        if outcome is not None:
            return outcome

        output_filename = f"{dataname}_TOM_{beta}.csv"
        produced_path = workdir / output_filename
        if not produced_path.exists():
            return TestOutcome(test_name, TestStatus.FAILED, f"missing output {output_filename}")

        produced = read_numeric_matrix(produced_path)
        reference = read_numeric_matrix(context.tests_dir / output_filename)
        is_equal, message = compare_matrices(produced, reference, TOM_TOLERANCE)

        status = TestStatus.PASSED if is_equal else TestStatus.FAILED
        return TestOutcome(test_name, status, message)


def run_build_graph_pow_tests(context, dataset_params):
    for dataset, gene_count, beta in dataset_params:
        yield timed(check_build_graph_pow, context, dataset, gene_count, beta)


def run_build_graph_tom_tests(context, dataset_params):
    for dataset, gene_count, beta in dataset_params:
        yield timed(check_build_graph_tom, context, dataset, gene_count, beta)


def count_build_graph_tests(dataset_params):
    return len(dataset_params)


# --- Palier 3: find_modules ----------------------------------------------------

def build_find_modules_command(tom_name, option, tau, l_value):
    """Return (program_args, produced_filename, reference_filename, test_name).

    produced_filename follows the subject literally (no "_groups_"); the
    downloaded reference files in data/tests/ still use the old
    "_groups_..." name (a bug in the reference generator, acknowledged by
    the teacher on Discord on 2026-07-14) - same content, different name.
    """
    if option == "static":
        program_args = [tom_name, "static", str(tau), str(l_value)]
        produced_filename = f"{tom_name}_static_{l_value}_{tau}.csv"
        reference_filename = f"{tom_name}_groups_static_{l_value}_{tau}.csv"
        test_name = f"find_modules {tom_name} static l={l_value} tau={tau}"
    else:
        program_args = [tom_name, option, str(tau)]
        produced_filename = f"{tom_name}_{option}_{tau}.csv"
        reference_filename = f"{tom_name}_groups_{option}_{tau}.csv"
        test_name = f"find_modules {tom_name} {option} tau={tau}"

    return program_args, produced_filename, reference_filename, test_name


def check_find_modules(context, tom_name, option, tau, l_value):
    program_args, produced_filename, reference_filename, test_name = build_find_modules_command(
        tom_name, option, tau, l_value,
    )

    with staged_directory(context.tests_dir, [tom_name]) as workdir:
        args = [str(workdir)] + program_args
        result, outcome = run_and_check(context, "find_modules.py", args, test_name)
        if outcome is not None:
            return outcome

        produced_path = workdir / produced_filename
        if not produced_path.exists():
            return TestOutcome(test_name, TestStatus.FAILED, f"missing output {produced_filename}")

        produced = read_groups_csv(produced_path)
        reference = read_groups_csv(context.tests_dir / reference_filename)
        is_equal, message = compare_groups(produced, reference)

        status = TestStatus.PASSED if is_equal else TestStatus.FAILED
        return TestOutcome(test_name, status, message)


def run_find_modules_tests(context, dataset_params):
    for dataset, gene_count, beta in dataset_params:
        tom_name = f"mrna2_{dataset}.csv_cleaned_{gene_count}.csv_TOM_{beta}.csv"
        for tau in TAU_VALUES:
            for l_value in STATIC_L_VALUES:
                yield timed(check_find_modules, context, tom_name, "static", tau, l_value)
            for option in ADAPTIVE_AND_DYNAMIC_OPTIONS:
                yield timed(check_find_modules, context, tom_name, option, tau, None)


def count_find_modules_tests(dataset_params):
    combinations_per_tom = len(STATIC_L_VALUES) + len(ADAPTIVE_AND_DYNAMIC_OPTIONS)
    return len(dataset_params) * len(TAU_VALUES) * combinations_per_tom


# --- Palier 4: extract_eigen ----------------------------------------------------

def check_extract_eigen(context, cleaned_name, group_name):
    test_name = f"extract_eigen {group_name}"

    with staged_directory(context.tests_dir, [cleaned_name, group_name]) as workdir:
        args = [str(workdir), cleaned_name, group_name]
        result, outcome = run_and_check(context, "extract_eigen.py", args, test_name)
        if outcome is not None:
            return outcome

        output_filename = f"{group_name}_eigen.csv"
        produced_path = workdir / output_filename
        if not produced_path.exists():
            return TestOutcome(test_name, TestStatus.FAILED, f"missing output {output_filename}")

        produced = read_eigen_csv(produced_path)
        reference = read_eigen_csv(context.tests_dir / output_filename)
        is_equal, message = compare_eigen(produced, reference, EIGEN_TOLERANCE)

        status = TestStatus.PASSED if is_equal else TestStatus.FAILED
        return TestOutcome(test_name, status, message)


def run_extract_eigen_tests(context, dataset_params):
    for dataset, gene_count, beta in dataset_params:
        cleaned_name = f"mrna2_{dataset}.csv_cleaned_{gene_count}.csv"
        tom_name = f"{cleaned_name}_TOM_{beta}.csv"
        for tau in TAU_VALUES:
            for l_value in STATIC_L_VALUES:
                group_name = f"{tom_name}_groups_static_{l_value}_{tau}.csv"
                yield timed(check_extract_eigen, context, cleaned_name, group_name)
            for option in ADAPTIVE_AND_DYNAMIC_OPTIONS:
                group_name = f"{tom_name}_groups_{option}_{tau}.csv"
                yield timed(check_extract_eigen, context, cleaned_name, group_name)


def count_extract_eigen_tests(dataset_params):
    combinations_per_tom = len(STATIC_L_VALUES) + len(ADAPTIVE_AND_DYNAMIC_OPTIONS)
    return len(dataset_params) * len(TAU_VALUES) * combinations_per_tom


# --- Reporting and CLI ------------------------------------------------------------

RUNNERS = [
    (1, "Palier 1 - clean_dataset", run_clean_dataset_tests, count_clean_dataset_tests),
    (2, "Palier 2 - build_graph POW", run_build_graph_pow_tests, count_build_graph_tests),
    (2, "Palier 2 - build_graph TOM", run_build_graph_tom_tests, count_build_graph_tests),
    (3, "Palier 3 - find_modules", run_find_modules_tests, count_find_modules_tests),
    (4, "Palier 4 - extract_eigen", run_extract_eigen_tests, count_extract_eigen_tests),
]


def filtered_dataset_params(datasets, gene_counts):
    return [
        (dataset, gene_count, beta)
        for (dataset, gene_count), beta in EXPECTED_BETA.items()
        if dataset in datasets and gene_count in gene_counts
    ]


def colorize(text, color, use_color):
    if not use_color:
        return text

    return f"{color}{text}{COLOR_RESET}"


def format_status_label(status, use_color):
    label = colorize(f"[{status.value}]", STATUS_COLORS[status.value], use_color)
    padding = " " * (STATUS_LABEL_WIDTH - len(status.value))
    return label + padding


def format_counts(counts, use_color):
    parts = []
    for status, count in counts.items():
        text = f"{status.value}={count}"
        color = STATUS_COLORS[status.value] if count else COLOR_GRAY
        parts.append(colorize(text, color, use_color))

    return " ".join(parts)


def format_duration(seconds):
    if seconds < 60:
        return f"{seconds:.1f}s"

    minutes, remaining_seconds = divmod(seconds, 60)
    return f"{int(minutes)}m{remaining_seconds:04.1f}s"


def truncate_text(text, max_length):
    if max_length <= 0:
        return ""
    if len(text) <= max_length:
        return text

    visible_length = max_length - len(TRUNCATION_SUFFIX)
    return text[:visible_length] + TRUNCATION_SUFFIX


def render_progress_bar(current, total, elapsed, last_outcome):
    """Render a single-line progress bar with count, ETA and the last result."""
    fraction = current / total if total else 1.0
    filled_width = int(PROGRESS_BAR_WIDTH * fraction)
    bar = "#" * filled_width + "-" * (PROGRESS_BAR_WIDTH - filled_width)
    percent = fraction * 100

    average_duration = elapsed / current if current else 0.0
    remaining_count = total - current
    eta_seconds = average_duration * remaining_count

    counters = f"{current}/{total} ({percent:5.1f}%)"
    timing = f"elapsed {format_duration(elapsed)} eta {format_duration(eta_seconds)}"
    prefix = f"[{bar}] {counters} {timing} "

    last_text = f"last: {last_outcome.status.value} {last_outcome.name}" if last_outcome else ""
    terminal_size = shutil.get_terminal_size(fallback=(FALLBACK_TERMINAL_WIDTH, 24))
    available_width = terminal_size.columns - len(prefix)
    last_text = truncate_text(last_text, available_width)

    return prefix + last_text


def print_test_line(outcome, position, total, context):
    status_label = format_status_label(outcome.status, context.use_color)
    duration_text = format_duration(outcome.duration)
    print(f"{status_label} [{position}/{total}] {outcome.name}: {outcome.message} ({duration_text})")


def print_outcomes(label, runner, counter, context, dataset_params):
    title = colorize(f"=== {label} ===", COLOR_BOLD + COLOR_CYAN, context.use_color)
    print(f"\n{title}")

    total = counter(dataset_params)
    show_progress_bar = context.is_tty and not context.verbose
    counts = {status: 0 for status in TestStatus}

    palier_start = time.perf_counter()
    last_outcome = None
    for position, outcome in enumerate(runner(context, dataset_params), start=1):
        counts[outcome.status] += 1
        last_outcome = outcome

        if context.verbose or outcome.status != TestStatus.PASSED:
            if show_progress_bar:
                print(f"\r{CLEAR_LINE}", end="")
            print_test_line(outcome, position, total, context)

        if show_progress_bar:
            elapsed = time.perf_counter() - palier_start
            bar_line = render_progress_bar(position, total, elapsed, last_outcome)
            print(f"\r{CLEAR_LINE}{bar_line}", end="", flush=True)

    if show_progress_bar:
        print(f"\r{CLEAR_LINE}", end="")

    palier_duration = time.perf_counter() - palier_start
    summary = format_counts(counts, context.use_color)
    print(f"--- {label}: {summary} ({format_duration(palier_duration)}) ---")

    return counts


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--src", type=Path, default=DEFAULT_SRC_DIR, help="directory with the .py programs")
    parser.add_argument("--tests", type=Path, default=DEFAULT_TESTS_DIR, help="directory with the reference files")
    parser.add_argument("--palier", type=int, choices=[1, 2, 3, 4], help="only run this palier")
    parser.add_argument("--datasets", nargs="+", choices=DATASETS, default=DATASETS)
    parser.add_argument("--gene-counts", nargs="+", type=int, choices=GENE_COUNTS, default=GENE_COUNTS)
    parser.add_argument("--quick", action="store_true", help="shortcut for --datasets small --gene-counts 50")
    parser.add_argument("--verbose", action="store_true", help="print PASS results too")
    parser.add_argument("--no-color", action="store_true", help="disable colored output")

    return parser.parse_args()


def main():
    args = parse_args()
    if args.quick:
        args.datasets = ["small"]
        args.gene_counts = [50]

    is_tty = sys.stdout.isatty()
    use_color = is_tty and not args.no_color
    context = TestContext(
        src_dir=args.src, tests_dir=args.tests, timeout=SUBPROCESS_TIMEOUT_SECONDS,
        verbose=args.verbose, use_color=use_color, is_tty=is_tty,
    )
    dataset_params = filtered_dataset_params(args.datasets, args.gene_counts)

    total_counts = {status: 0 for status in TestStatus}
    run_start = time.perf_counter()
    for palier, label, runner, counter in RUNNERS:
        if args.palier is not None and palier != args.palier:
            continue

        counts = print_outcomes(label, runner, counter, context, dataset_params)
        for status, count in counts.items():
            total_counts[status] += count

    run_duration = time.perf_counter() - run_start
    summary = format_counts(total_counts, use_color)
    separator = colorize("=" * 60, COLOR_BOLD, use_color)
    print(f"\n{separator}")
    print(colorize("TOTAL", COLOR_BOLD, use_color) + f": {summary} ({format_duration(run_duration)})")
    print(separator)

    if total_counts[TestStatus.FAILED] or total_counts[TestStatus.ERROR]:
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
