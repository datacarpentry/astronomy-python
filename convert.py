import nbformat as nbf
from glob import glob
#from pypandoc import convert_text
import re
from textwrap import wrap

# Collect a list of all notebooks in the content folder
notebooks = glob("*.ipynb")

# Text to look for in adding tags
tag_search_dict = {
    "remove-cell": "#hide\n",
    "hide-cell": "#hide\n",
    "remove-input": "#hide_input\n",
    "hide-input": "#hide_input\n",
    "remove-output": "#hide_output\n",
    "hide-output": "#hide_output\n",
}

WRAP_OPTIONS = dict(break_long_words=False,
                    break_on_hyphens=False)

def wrap_source(source):
    res = []
    for line in source.split('\n'):
        if len(line) > 80:
            line = '\n'.join(wrap(line, **WRAP_OPTIONS))
        res.append(line)
    return '\n'.join(res)

def remove_first_line(text):
    _, _, rest = text.partition('\n')
    return rest

def indent_text(text, prefix):
    res = []
    for line in text.split('\n'):
        res.append(prefix + line + '\n')
    return ''.join(res)

def format_exercise(cell):
    res = []
    res.append('> ## Exercise\n')
    text = remove_first_line(cell['source'])
    res.append(indent_text(text, '> '))
    return ''.join(res)


def process_markdown(cell):
    if cell['source'].startswith('### Exercise'):
        cell['source'] = format_exercise(cell)
    else:
        cell['source'] = re.sub('\n>', '\n> :', cell['source'])
    cell['source'] = wrap_source(cell['source'])


def process_code(cell):
    """Translate code blocks
    """
    text = cell['source']

    # If it's a solution, remove the first line
    is_solution = text.startswith('# Solution')
    if is_solution:
        text = remove_first_line(text)

    # Assemble the exercise formatting
    cell['cell_type'] = 'raw'
    res = []
    res.append('\n~~~')
    res.append(text)
    res.append('~~~')
    res.append('{: .language-python}')

    if cell['outputs']:
        res.append('\n~~~')
        for output in cell['outputs']:
            try:
                res.append(output['text'])
            except KeyError:
                """Ignoring non-text output"""
                pass
        res.append('~~~')
        res.append('{: .output}\n')

    cell['source'] = '\n'.join(res)

    # If it's a solution, indent it and add the suffix
    if is_solution:
        cell['source'] = indent_text(cell['source'], '> > ')
        cell['source'] += '> {: .solution}\n'
        cell['source'] += '{: .challenge}\n'

def process_cell(cell):
    if cell['cell_type'] == 'raw':
        return

    cell_tags = cell.get('metadata', {}).get('tags', [])
    if 'remove-cell' in cell_tags:
        cell['cell_type'] = 'raw'
        cell['source'] = ''
        cell['outputs'] = []
        return

    if cell['cell_type'] == 'markdown':
        process_markdown(cell)

    if cell['cell_type'] == 'code':
        process_code(cell)

def process_notebook(path):
    """
    """
    ntbk = nbf.read(path, nbf.NO_CONVERT)

    for cell in ntbk.cells:
        process_cell(cell)

    nbf.write(ntbk, path)

# Search through each notebook and look for the text, add a tag if necessary
for path in notebooks:
    print('converting', path)
    process_notebook(path)
