"""
wrap_project_tex.py — Standardizes LaTeX Line Wrapping for Paper IV
Format Utility: Strict 80-Character Width Standardizer with Math-Preservation

This utility formats the entire LaTeX manuscript suite to a standard 80-character 
line wrapping width. It optimizes the source text layout for human readability and 
version control tracking while ensuring that code blocks, math environments, and tables
are strictly preserved.
"""
import os
import textwrap

def split_latex_line(line):
    """Preprocess a raw LaTeX line, splitting inline subheadings and list items.

    This function prevents formatting issues where content following '\\paragraph'
    or '\\item' commands is joined onto the same line as the command.

    Args:
        line: The raw input string line from the LaTeX file.

    Returns:
        A list of strings containing one or two lines depending on command matches.
    """
    stripped = line.strip()
    # Handle \paragraph{Title} Text...
    if stripped.startswith('\\paragraph{'):
        brace_count = 0
        idx = -1
        for i, char in enumerate(stripped):
            if char == '{':
                brace_count += 1
            elif char == '}':
                brace_count -= 1
                if brace_count == 0:
                    idx = i
                    break
        if idx != -1 and idx + 1 < len(stripped):
            heading = stripped[:idx+1]
            rest = stripped[idx+1:].strip()
            if rest:
                # Retain original indentation of the line for the heading
                indent = line[:line.find('\\paragraph{')]
                return [indent + heading, indent + rest]
                
    # Handle \item Text...
    elif stripped.startswith('\\item '):
        # Retain original indentation of the line for \item
        indent = line[:line.find('\\item ')]
        rest = stripped[5:].strip()
        if rest:
            return [indent + '\\item', indent + rest]
            
    return [line]

def wrap_paragraph(paragraph_lines, width=80):
    """Enforce a strict 80-character textwrap on a list of accumulated paragraph lines.

    Args:
        paragraph_lines: List of accumulated raw string lines forming a paragraph.
        width: Maximum target character width for the wrapping (default: 80).

    Returns:
        List of formatted lines wrapped within the width constraint.
    """
    if not paragraph_lines:
        return []
    joined_text = " ".join(line.strip() for line in paragraph_lines)
    wrapped_lines = textwrap.wrap(joined_text, width=width, break_long_words=False, break_on_hyphens=False)
    return wrapped_lines

def process_file_content(content):
    """Standardize the line wrapping across the text content of a LaTeX file.

    This engine parses the document linearly, buffering normal paragraph lines
    and wrapping them to 80 characters, while ignoring specific environments
    (e.g., equations, aligns, tabulars, figures, tables) and single-line commands
    to prevent syntax breakage.

    Args:
        content: The raw text string of the LaTeX file.

    Returns:
        The fully formatted text string with CRLF newlines.
    """
    # Preprocess lines to split headings and items
    raw_lines = content.splitlines()
    lines = []
    for rl in raw_lines:
        lines.extend(split_latex_line(rl))
        
    output_lines = []
    current_paragraph = []
    
    NO_WRAP_ENVIRONMENTS = {'equation', 'equation*', 'align', 'align*', 'gather', 'gather*', 'tabular', 'tabularx', 'table', 'figure'}
    
    DONT_WRAP_COMMANDS = (
        '\\section', '\\subsection', '\\subsubsection', '\\paragraph',
        '\\label', '\\caption', '\\input', '\\usepackage', '\\documentclass',
        '\\bibliographystyle', '\\bibliography', '\\maketitle', '\\author',
        '\\title', '\\date', '\\item', '\\newpage', '\\clearpage', '\\appendix',
        '\\centering', '\\includegraphics', '\\hline', '\\toprule', '\\midrule',
        '\\bottomrule', '\\newcolumntype', '\\newtheorem', '\\def', '\\newcommand',
        '\\renewcommand', '\\footnote', '\\begin', '\\end'
    )
    
    active_no_wrap_envs = []
    
    for line in lines:
        stripped = line.strip()
        
        is_begin = stripped.startswith('\\begin{')
        is_end = stripped.startswith('\\end{')
        
        env_name = ""
        if is_begin:
            try:
                env_name = stripped.split('{')[1].split('}')[0]
            except:
                pass
        elif is_end:
            try:
                env_name = stripped.split('{')[1].split('}')[0]
            except:
                pass
        
        if is_begin and env_name in NO_WRAP_ENVIRONMENTS:
            active_no_wrap_envs.append(env_name)
        
        in_no_wrap = len(active_no_wrap_envs) > 0
        
        should_not_wrap_line = (
            stripped == "" or
            stripped.startswith('%') or
            stripped.startswith('\\begin') or
            stripped.startswith('\\end') or
            stripped.startswith(DONT_WRAP_COMMANDS) or
            in_no_wrap
        )
        
        if should_not_wrap_line:
            if current_paragraph:
                output_lines.extend(wrap_paragraph(current_paragraph))
                current_paragraph = []
            output_lines.append(line)
        else:
            # Handle lines ending with \\ or containing \newline by flushing them immediately
            if stripped.endswith('\\\\') or '\\newline' in stripped:
                current_paragraph.append(line)
                output_lines.extend(wrap_paragraph(current_paragraph))
                current_paragraph = []
            else:
                current_paragraph.append(line)
            
        if is_end and env_name in NO_WRAP_ENVIRONMENTS:
            if active_no_wrap_envs and active_no_wrap_envs[-1] == env_name:
                active_no_wrap_envs.pop()
                
    if current_paragraph:
        output_lines.extend(wrap_paragraph(current_paragraph))
        
    return "\r\n".join(output_lines) + "\r\n"

def main():
    """Main execution entry point to gather and format all project LaTeX files."""
    base_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.join(base_dir, "..", "..")
    sections_dir = os.path.join(project_dir, "manuscript", "sections")
    
    files_to_format = []
    
    # Gather all .tex files in manuscript/sections
    if os.path.exists(sections_dir):
        for file in os.listdir(sections_dir):
            if file.endswith(".tex"):
                files_to_format.append(os.path.join(sections_dir, file))
                
    # Gather main.tex and supplemental.tex
    main_tex = os.path.join(project_dir, "manuscript", "main.tex")
    if os.path.exists(main_tex):
        files_to_format.append(main_tex)
        
    supp_tex = os.path.join(project_dir, "supplements", "supplemental.tex")
    if os.path.exists(supp_tex):
        files_to_format.append(supp_tex)
        
    print(f"Found {len(files_to_format)} files to format.")
    
    for path in files_to_format:
        print(f"Formatting {os.path.basename(path)}...")
        with open(path, "rb") as f:
            content = f.read().decode("utf-8")
        
        formatted = process_file_content(content)
        
        with open(path, "wb") as f:
            f.write(formatted.encode("utf-8"))
            
    print("Global project formatting completed successfully!")

if __name__ == "__main__":
    main()
