#!/usr/bin/env python3
"""
check_orphaned_macros.py — LaTeX Macro Usage Analyzer for Paper IV
Verification Utility: Dynamic Variables Orphan Detection and Usage Statistics

This utility analyzes the complete LaTeX manuscript suite (manuscript + supplements)
to identify orphaned macro definitions in dynamic_variables.tex. It provides
comprehensive usage statistics for all defined macros, helping maintain clean
and efficient variable files by flagging unused definitions that can be safely
removed.

The script scans all .tex files recursively, using regex pattern matching to detect
macro invocations while excluding the definition file itself to avoid false positives.

Usage:
    python scripts/utility/check_orphaned_macros.py

Output:
    - Total macro count from dynamic_variables.tex
    - List of completely orphaned macros (usage count = 0)
    - List of rarely-used macros (usage count = 1)
    - Statistics summary (orphaned, single-use, multi-use counts)

Author: Generated during Paper IV audit (2026-05-27)
"""

import os
import re


def extract_defined_macros(dyn_file_path):
    """Extract all macro names defined in dynamic_variables.tex.

    Args:
        dyn_file_path: Path to the dynamic_variables.tex file.

    Returns:
        List of macro names (without backslash prefix) defined via \\newcommand.
    """
    with open(dyn_file_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # Pattern: \newcommand{\MacroName}{value}
    macro_names = re.findall(r'\\newcommand\{\\([a-zA-Z]+)\}', content)
    return macro_names


def collect_tex_files(root_dir, exclude_filename='dynamic_variables.tex'):
    """Recursively collect all .tex files except the definition file.

    Args:
        root_dir: Root directory to start the recursive search.
        exclude_filename: Filename to exclude (default: dynamic_variables.tex).

    Returns:
        List of full paths to .tex files found, excluding the definition file.
    """
    tex_files = []
    for root, dirs, files in os.walk(root_dir):
        for f in files:
            if f.endswith('.tex') and f != exclude_filename:
                tex_files.append(os.path.join(root, f))
    return tex_files


def count_macro_usage(macro_names, tex_files):
    """Count how many times each macro is used across all .tex files.

    Args:
        macro_names: List of macro names to search for.
        tex_files: List of .tex file paths to scan.

    Returns:
        Dictionary mapping macro names to usage counts {macro_name: count}.
    """
    usage_counts = {}

    for macro in macro_names:
        count = 0
        for tex_file in tex_files:
            try:
                with open(tex_file, 'r', encoding='utf-8') as f:
                    content = f.read()

                # Pattern: \MacroName (not followed by letter to avoid partial matches)
                pattern = r'\\' + re.escape(macro) + r'(?![a-zA-Z])'
                matches = re.findall(pattern, content)
                count += len(matches)

            except Exception as e:
                print(f'Error reading {tex_file}: {e}')

        usage_counts[macro] = count

    return usage_counts


def categorize_macros(usage_counts):
    """Categorize macros by usage frequency.

    Args:
        usage_counts: Dictionary of {macro_name: usage_count}.

    Returns:
        Tuple of three sorted lists: (orphaned, used_once, used_multiple).
        - orphaned: macros with count == 0
        - used_once: macros with count == 1
        - used_multiple: macros with count > 1
    """
    orphaned = sorted([m for m, c in usage_counts.items() if c == 0])
    used_once = sorted([m for m, c in usage_counts.items() if c == 1])
    used_multiple = sorted([m for m, c in usage_counts.items() if c > 1])

    return orphaned, used_once, used_multiple


def print_results(orphaned, used_once, used_multiple, usage_counts):
    """Print formatted analysis results to stdout.

    Args:
        orphaned: List of orphaned macro names.
        used_once: List of single-use macro names.
        used_multiple: List of multi-use macro names.
        usage_counts: Dictionary of {macro_name: usage_count}.
    """
    print(f'=== RESULTS ===')
    print(f'Orphaned macros (never used): {len(orphaned)}')
    print(f'Macros used once: {len(used_once)}')
    print(f'Macros used multiple times: {len(used_multiple)}\n')

    if orphaned:
        print('COMPLETELY ORPHANED MACROS:')
        for macro in orphaned:
            print(f'  \\{macro}')
        print()

    if used_once:
        print('MACROS USED ONLY ONCE (verify if necessary):')
        display_limit = 10
        for macro in used_once[:display_limit]:
            print(f'  \\{macro} (used {usage_counts[macro]} time)')

        if len(used_once) > display_limit:
            remaining = len(used_once) - display_limit
            print(f'  ... and {remaining} more macros')


def main():
    """Main execution: analyze dynamic_variables.tex macro usage."""
    # File paths relative to project root
    dyn_file = r'manuscript\dynamic_variables.tex'
    project_root = '.'

    # Step 1: Extract defined macros
    defined_macros = extract_defined_macros(dyn_file)
    print(f'Found {len(defined_macros)} macros defined in dynamic_variables.tex\n')

    # Step 2: Collect all .tex files (excluding definition file)
    tex_files = collect_tex_files(project_root, exclude_filename='dynamic_variables.tex')
    print(f'Scanning {len(tex_files)} .tex files')
    print(f'Files: {[os.path.basename(f) for f in tex_files[:5]]}...\n')

    # Step 3: Count usage across all files
    usage_counts = count_macro_usage(defined_macros, tex_files)

    # Step 4: Categorize by frequency
    orphaned, used_once, used_multiple = categorize_macros(usage_counts)

    # Step 5: Print results
    print_results(orphaned, used_once, used_multiple, usage_counts)


if __name__ == '__main__':
    main()
