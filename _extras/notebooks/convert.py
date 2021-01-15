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

def wrap_text(text):
    """Wrap lines so they don't exceed 80 chars.

    text: string of newline-separated lines

    returns: string of newline-separated lines
    """
    res = []
    for line in text.split('\n'):
        if len(line) > 80:
            wrapped = wrap(line, **WRAP_OPTIONS)
            res.extend(wrapped)
        else:
            res.append(line)
    return '\n'.join(res)

def remove_first_line(text):
    """Remove the first line from a text block.

    source: string of newline-separated lines

    returns: string of newline-separated lines
    """
    _, _, rest = text.partition('\n')
    return rest

def indent_text(text, prefix):
    res = []
    for line in text.split('\n'):
        res.append(prefix + line)
    return '\n'.join(res)

def format_exercise(cell):
    res = []
    res.append('> ## Exercise\n')
    text = remove_first_line(cell['source'])
    text = wrap_text(text)
    #print('--------------')
    #print(text)
    #print('--------------')
    text = indent_text(text, '> ')
    #print('--------------')
    #print(text)
    #print('--------------')
    res.append(text)
    joined = ''.join(res)
    #print('--------------')
    #print(joined)
    #print('--------------')
    return joined

def process_markdown(cell):
    if cell['source'].startswith('### Exercise'):
        cell['source'] = format_exercise(cell)
    else:
        cell['source'] = re.sub('\n>', '\n> :', cell['source'])
        cell['source'] = wrap_text(cell['source'])

def truncate_text(text, cutoff=11):
    """Cut off outputs that exceed 11 lines.
    """
    t = text.split('\n')
    if len(t) > cutoff:
        t = t[:cutoff]
        t.append('[Output truncated]')
        return '\n'.join(t)
    else:
        return text

def add_text(output, res):
    try:
        text = output['text']
        res.append(truncate_text(text))
    except KeyError:
        """Ignoring non-text output"""
        pass

def add_data(output, res):
    try:
        data = output['data']
        text = data['text/plain']
        res.append(truncate_text(text))
    except KeyError:
        """Ignoring non-data output"""
        pass

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

    if cell['outputs'] and not is_solution:
        res.append('\n~~~')
        for output in cell['outputs']:
            add_text(output, res)
            add_data(output, res)

        res.append('~~~')
        res.append('{: .output}\n')

    cell['source'] = '\n'.join(res)

    # If it's a solution, indent it and add the suffix
    if is_solution:
        text = indent_text(cell['source'], '> > ')
        res = []
        res.append('>')
        res.append('> > ## Solution')
        res.append(text)
        res.append('> {: .solution}')
        res.append('{: .challenge}')
        res.append('\n')

    cell['source'] = '\n'.join(res)

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
