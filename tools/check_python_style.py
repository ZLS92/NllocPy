# =============================================================================
# Author: Luigi Sante Zampa
# Affiliation: NllocPy project (institutional affiliation not specified)
# Created: 2026-11-04
# Last modified: 2026-11-10
#
# Purpose:
#   Audit Python sources against NllocPy project conventions. The core modules
#   are enforced by CI; legacy modules can be audited in advisory mode while
#   behavioral regression tests are added before further refactoring.
# =============================================================================

# -----------------------------------------------------------------------------
# Import external libraries
# -----------------------------------------------------------------------------
import argparse
import ast
from pathlib import Path


# -----------------------------------------------------------------------------
# Shared audit configuration
# -----------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]
STRICT_FILES = (
    "NllocPy/grid_tools.py",
    "NllocPy/model_tools.py",
    "NllocPy/table_tools.py",
)
HEADER_FIELDS = (
    "# Author:",
    "# Affiliation:",
    "# Created:",
    "# Last modified:",
    "# Purpose:",
)


# =============================================================================
# Standalone helpers for Python style and structure inspection
# =============================================================================

# -----------------------------------------------------------------------------
def audit_module(path, *, strict=False):
    """Collect code-style problems in an individual Python source file.

    Parameters
    ----------
    path : str or pathlib.Path
        Source file to inspect.
    strict : bool, default False
        Enable all structural, documentation and nested-function checks.

    Returns
    -------
    list of str
        Human-readable compliance findings. An empty list means no findings.

    Examples
    --------
    >>> isinstance(audit_module(ROOT / "NllocPy/table_tools.py"), list)
    True
    """
    path = Path(path)
    content = path.read_text(encoding="utf-8")
    lines = content.splitlines()
    problems = []

    # Basic provenance and version-header checks are useful on every module.
    first_lines = "\n".join(lines[:30])

    for field in HEADER_FIELDS:
        if field not in first_lines:
            problems.append(f"Missing header field {field}")

    try:
        tree = ast.parse(content, filename=str(path))
    except SyntaxError as exc:
        problems.append(f"Invalid Python syntax: {exc}")
        return problems

    if not strict:
        output = problems
        return output

    for node in tree.body:
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            continue

        # Each standalone function/class must have a nearby separator.
        preceding = "\n".join(lines[max(0, node.lineno - 5):node.lineno - 1])

        if "# ---" not in preceding:
            problems.append(f"L{node.lineno}: missing separator before {node.name}")

        doc = ast.get_docstring(node)

        if not doc:
            problems.append(f"L{node.lineno}: missing docstring for {node.name}")
        else:
            # Classes describe their constructor parameters; they do not
            # return values, so a Returns section is not appropriate.
            required = ("Parameters", "Examples")
            if not isinstance(node, ast.ClassDef):
                required = ("Parameters", "Returns", "Examples")

            if not all(part in doc for part in required):
                problems.append(f"L{node.lineno}: incomplete docstring for {node.name}")

        for inner in ast.walk(node):
            if inner is node:
                continue

            if isinstance(inner, (ast.FunctionDef, ast.AsyncFunctionDef)):
                if isinstance(node, ast.ClassDef) and inner.name == "__init__":
                    continue
                if isinstance(node, ast.ClassDef):
                    problems.append(f"L{inner.lineno}: class method {inner.name} must be external")
                else:
                    problems.append(f"L{inner.lineno}: nested function {inner.name}")

            if isinstance(inner, (ast.Import, ast.ImportFrom)):
                problems.append(f"L{inner.lineno}: import must be at module header")

            if isinstance(inner, ast.Return) and isinstance(inner.value, ast.Tuple):
                problems.append(f"L{inner.lineno}: assign multiple return values before returning")

    output = problems
    return output


# -----------------------------------------------------------------------------
def audit_paths(paths, *, strict=False):
    """Check multiple Python sources and collect their findings.

    Parameters
    ----------
    paths : sequence of str or pathlib.Path
        Files to analyze.
    strict : bool
        Whether to apply the complete structural checks.

    Returns
    -------
    dict
        Dictionary mapping each file path to a list of findings.

    Examples
    --------
    >>> isinstance(audit_paths([ROOT / "NllocPy/grid_tools.py"]), dict)
    True
    """
    findings = {}

    for path in paths:
        problems = audit_module(path, strict=strict)
        findings[str(path)] = problems

    output = findings
    return output


# =============================================================================
# Command-line interface
# =============================================================================

# -----------------------------------------------------------------------------
def main():
    """Run the strict or advisory convention audit from the command line.

    Parameters
    ----------
    None
        Arguments are parsed from the Python command line.

    Returns
    -------
    int
        0 when strict selected files pass, or when using advisory --all;
        1 when strict selected files have findings.

    Examples
    --------
    >>> # python tools/check_python_style.py --strict
    >>> # python tools/check_python_style.py --all
    """
    parser = argparse.ArgumentParser(
        description="Check NllocPy Python code conventions"
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Check the modernized modules with all structural rules",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Advisory report of full conventions across legacy source modules",
    )

    args = parser.parse_args()

    if args.all:
        paths = sorted((ROOT / "NllocPy").glob("*.py"))
        paths.append(ROOT / "setup.py")
        strict = True  # Full audit, advisory exit status for legacy migration.
    else:
        paths = [ROOT / path for path in STRICT_FILES]
        strict = True

    results = audit_paths(paths, strict=strict)
    errors = 0

    for path, problems in results.items():
        name = str(Path(path).relative_to(ROOT))

        if problems:
            for problem in problems:
                print(f"{name}: {problem}")
                errors += 1
        else:
            print(f"{name}: OK")

    if args.all:
        print(f"Advisory full-package findings: {errors}")
        return 0

    output = int(errors > 0)
    return output


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    raise SystemExit(main())
