import os
import zipfile
import re

base_dir = r'c:\Users\ricca\Documenti\Progetti IA\branching papers\paper4-incommensurability'
out_dir = r'c:\Users\ricca\Documenti\Progetti IA\branching papers\Arxiv\paper4'
os.makedirs(out_dir, exist_ok=True)
zip_path = os.path.join(out_dir, 'paper4_arxiv_full.zip')

def resolve_tex_content(filepath, replacements):
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    for old, new in replacements.items():
        content = content.replace(old, new)
    return content

def add_folder(z, folder_path, prefix='', allowed_exts=None):
    for root, dirs, files in os.walk(folder_path):
        if '__pycache__' in root: continue
        for file in files:
            ext = os.path.splitext(file)[1].lower()
            if allowed_exts is None or ext in allowed_exts:
                full_path = os.path.join(root, file)
                arc_name = os.path.join(prefix, os.path.relpath(full_path, folder_path))
                arc_name = arc_name.replace('\\', '/')
                z.write(full_path, arc_name)

with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
    
    # 1. Main Tex
    main_path = os.path.join(base_dir, 'manuscript', 'main.tex')
    
    # 2. Add Pre-compiled Supplemental PDF to anc/
    supp_pdf_path = os.path.join(base_dir, 'output', 'Supplemental Material - The Incommensurability Principle.pdf')
    if os.path.exists(supp_pdf_path):
        z.write(supp_pdf_path, 'anc/Supplemental_Material.pdf')
    else:
        print("Warning: Supplemental PDF not found!")
    
    # Remove \externaldocument from main.tex to avoid ?? errors as much as possible
    # We replace it with an empty string so we don't accidentally comment out \AtBeginDocument on the same line!
    main_rep = {
        r'\externaldocument[S-]{../supplements/supplemental}': r''
    }
    z.writestr('main.tex', resolve_tex_content(main_path, main_rep))
    
    # 3. Dynamic Variables & References
    z.write(os.path.join(base_dir, 'manuscript', 'dynamic_variables.tex'), 'dynamic_variables.tex')
    z.write(os.path.join(base_dir, 'manuscript', 'references.bib'), 'references.bib')
    
    # Optional .bbl files
    main_bbl = os.path.join(base_dir, 'manuscript', 'main.bbl')
    if os.path.exists(main_bbl):
        z.write(main_bbl, 'main.bbl')
        
    supp_bbl = os.path.join(base_dir, 'supplements', 'supplemental.bbl')
    if os.path.exists(supp_bbl):
        z.write(supp_bbl, 'supplemental.bbl')
    
    # 4. Folders
    tex_exts = ['.tex', '.dat', '.png', '.jpg', '.pdf', '.bib']
    add_folder(z, os.path.join(base_dir, 'manuscript', 'sections'), 'sections', tex_exts)
    add_folder(z, os.path.join(base_dir, 'manuscript', 'figures'), 'figures', tex_exts)
    
    # 5. Ancillary Files (Scripts & README)
    readme_path = os.path.join(base_dir, 'README.md')
    if os.path.exists(readme_path):
        z.write(readme_path, 'anc/README.md')
        
    add_folder(z, os.path.join(base_dir, 'scripts'), 'anc/scripts')

print(f"Archivio salvato con successo: {zip_path}")
