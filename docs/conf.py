import git

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'PyOrbb'
copyright = '2025, TheoCheM VU Amsterdam'
author = 'TheoCheM VU Amsterdam'

# get release information
repo = git.Repo('..')

tags = sorted(repo.tags, key=lambda t: t.commit.committed_datetime)
if len(tags) == 0:
    latest_tag = None
    release = 'vUnknown'
else:
    latest_tag = tags[-1]
    release = latest_tag.name

print('Git data:')
print('\tRepository:    ', repo)
print('\tHeads:         ', repo.heads)
print('\tTags:          ', tags)
print('\tLatest Tag:    ', repr(latest_tag))
print('\tLatest Version:', release)

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration',
    "sphinx.ext.autodoc",
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    # 'sphinx.ext.autosummary',
    "sphinx_autodoc_typehints",
    "sphinx_click",
    "sphinx_design",
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


autodoc_default_options = {
    'autosummary': False,
}

modindex_common_prefix = ['template.']

html_theme_options = {
  # "show_nav_level": 2,
  # "navigation_depth": 2,
  "navbar_end": ["star"],
  "navbar_center": [],
  "navbar_start": ['logo'],
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_favicon = 'https://avatars.githubusercontent.com/u/119413491' 
html_theme = 'pydata_sphinx_theme'  # pip install pydata-sphinx-theme
html_static_path = ['_static']
add_module_names = False

# custom variables
rst_epilog = f"""
.. |ProjectName| replace:: {project}
.. |ProjectVersion| replace:: {release}
.. |eV| replace:: :math:`\\text{{eV}}`
.. |kcal/mol| replace:: :math:`\\text{{kcal mol}}^{{-1}}`
.. |Orbitals| replace:: :class:`~pyfmo.orbitals.objects.Orbitals`
.. |Orbital| replace:: :class:`~pyfmo.orbitals.objects.Orbital`
.. |MOs| replace:: :class:`~pyfmo.orbitals.objects.MOs`
.. |MO| replace:: :class:`~pyfmo.orbitals.objects.MO`
.. |SFOs| replace:: :class:`~pyfmo.orbitals.objects.SFOs`
.. |SFO| replace:: :class:`~pyfmo.orbitals.objects.SFO`
"""
