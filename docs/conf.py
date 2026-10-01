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
    "sphinx_autodoc_typehints",
    "sphinxarg.ext",
    "sphinx_design",
    "sphinx_tabs.tabs",
    "sphinxcontrib.video",
]

hoverxref_auto_ref = True
hoverxref_role_types = {
    'ref': 'tooltip',   # or 'modal' for a popup box
}
# enable external link icons for "sphinx_new_tab_link":
new_tab_link_show_external_link_icon = True

# configuration of "sphinxcontrib.lightbox2"
lightbox2_image_fade_duration = 300
lightbox2_resize_duration = 300
lightbox2_disable_scrolling = True
home_page_in_toc = True
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


autodoc_default_options = {
    'autosummary': False,
}

modindex_common_prefix = ['template.']

html_theme_options = {
  # "show_toc_level": 2,
  "show_nav_level": 2,
  # "navigation_depth": 2,
  # "navbar_end": ["star"],
  "navbar_center": [],
  # "navbar_start": ['logo'],
"logo": {
    # In a left-to-right context, screen readers will read the alt text
    # first, then the text, so this example will be read as "P-G-G-P-Y
    # (short pause) Home A pretty good geometry package"
    "text": f"PyOrbb {release} documentation",
    "image_light": "_static/images/icon_12.png",
    "image_dark": "_static/images/icon_12.png",
    }
}

html_js_files = ['gui_tabs.js']

toc_object_entries_show_parents = 'all'
# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_favicon = 'https://avatars.githubusercontent.com/u/119413491' 
html_theme = 'pydata_sphinx_theme'  # pip install pydata-sphinx-theme
html_static_path = ['_static']
add_module_names = False
toc_object_entries_show_parents = 'all'
# html_sidebars = { '**': ['globaltoc.html'] }

# custom variables
rst_epilog = f"""
.. |ProjectName| replace:: {project} 
.. |logo| image:: /_static/images/icon_12.png
    :height: 70px
    :align: middle
    :class: no-scaled-link
.. |ProjectVersion| replace:: {release}
.. |eV| replace:: :math:`\\text{{eV}}`
.. |kcal/mol| replace:: :math:`\\text{{kcal mol}}^{{-1}}`
.. |Orbitals| replace:: :class:`~pyorbb.orbitals.objects.Orbitals`
.. |Orbital| replace:: :class:`~pyorbb.orbitals.objects.Orbital`
.. |MOs| replace:: :class:`~pyorbb.orbitals.objects.MOs`
.. |MO| replace:: :class:`~pyorbb.orbitals.objects.MO`
.. |FMOs| replace:: :class:`~pyorbb.orbitals.objects.FMOs`
.. |FMO| replace:: :class:`~pyorbb.orbitals.objects.FMO`
.. |PyOrbb| replace:: :program:`PyOrbb`
.. |AMS| replace:: :program:`AMS`
.. |ADF| replace:: :program:`ADF`
.. |densf| replace:: :program:`densf`
.. |main art| replace:: Yuman Hordijk, Steven E. Beutick, Xiaobo Sun, Laurens Groot, Tori Gijzen, Jordi Poater, Trevor A. Hamlin, F. Matthias Bickelhaupt, Célia Fonseca Guerra, "PyOrbb – Automated Analyses of Orbital-Interaction Mechanisms" *Journal of Computational Chemistry*, **2026**, e70512, DOI=https://doi.org/10.1002/jcc.70512.
"""
