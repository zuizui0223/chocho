"""Check generated figure inventory; distinguish raster identity from numerical reconstruction."""
import hashlib,json,platform
from pathlib import Path
import numpy,scipy,matplotlib
P=Path(__file__).resolve().parent
m=json.loads((P/'word_candidate_manifest.json').read_text())
rows=[]
for f in m['figures']:
    stem=f['renderer_stem']
    png=P/(stem+'.png'); pdf=P/(stem+'.pdf')
    b=png.read_bytes(); v=pdf.read_bytes()
    assert b.startswith(b'\x89PNG\r\n\x1a\n') and len(b)>1000
    assert v.startswith(b'%PDF-') and len(v)>1000
    h=hashlib.sha256(b).hexdigest()
    rows.append({'figure':f['number'],'renderer_stem':stem,'png_sha256':h,
                 'reference_png_match':h==f['embedded_png_sha256'],
                 'pdf_sha256':hashlib.sha256(v).hexdigest()})
assert len(rows)==8 and len({r['figure'] for r in rows})==8
report={'status':'EIGHT_FIGURES_RENDERED','python':platform.python_version(),
        'numpy':numpy.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,
        'exact_reference_raster_matches':sum(r['reference_png_match'] for r in rows),'figures':rows,
        'scope':'Reference PNG mismatch may reflect renderer environment; it is not by itself evidence of a numerical difference. Full raw-model reconstruction and reviewer release are separate gates.'}
(P/'render_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
