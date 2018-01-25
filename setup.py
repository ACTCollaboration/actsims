
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': '',
    'author': 'Alexander van Engelen',
    'version': 'beta',
    'install_requires': [],
    'packages': ['actsims'
                ],
    'scripts': [],
    'name': 'actsims'
}

setup(**config)
