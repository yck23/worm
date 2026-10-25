"""Execute the response notebook with NumPy/Matplotlib/IPython; save its outputs."""
from pathlib import Path
import io
import json
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from IPython.core.interactiveshell import InteractiveShell
from IPython.display import Image, display
from IPython.utils.capture import capture_output

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / 'results/response_measurements/response_measurements.ipynb'

def show_figures(*args, **kwargs):
    for number in plt.get_fignums():
        stream = io.BytesIO()
        plt.figure(number).savefig(stream, format='png', dpi=120)
        display(Image(data=stream.getvalue()))
    plt.close('all')


def main():
    notebook = json.loads(PATH.read_text(encoding='utf-8'))
    os.chdir(PATH.parent)
    InteractiveShell.instance()
    plt.show = show_figures
    namespace = {'__name__': '__main__'}
    execution = 0
    for cell in notebook['cells']:
        if cell['cell_type'] != 'code':
            continue
        execution += 1
        source = ''.join(cell['source'])
        if 'def ' in source or 'lambda ' in source:
            raise ValueError('Keep the notebook as direct analysis cells, without custom functions')
        with capture_output() as captured:
            exec(compile(source, str(PATH) + ':' + cell['id'], 'exec'), namespace)
        outputs = []
        if captured.stdout:
            outputs.append(dict(output_type='stream', name='stdout', text=captured.stdout.splitlines(keepends=True)))
        if captured.stderr:
            outputs.append(dict(output_type='stream', name='stderr', text=captured.stderr.splitlines(keepends=True)))
        for item in captured.outputs:
            outputs.append(dict(output_type='display_data', data=item.data, metadata=item.metadata))
        cell['execution_count'] = execution
        cell['outputs'] = outputs
        print('Executed', cell['id'], flush=True)
    notebook['metadata']['language_info']['version'] = sys.version.split()[0]
    PATH.write_text(json.dumps(notebook, indent=1) + '\n', encoding='utf-8')
    print(f'Saved {execution} executed code cells: {PATH}')

if __name__ == '__main__':
    main()
