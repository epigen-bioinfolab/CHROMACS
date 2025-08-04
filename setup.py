from setuptools import setup, find_packages

setup(
    name="chromacs",
    version="1.0",
    packages=find_packages(),
    install_requires=[
        'pandas',
        'pyyaml',
    ],
    include_package_data=True,
    package_data={
        'chromacs': [
            '*.R',
            'assets/*.png',
            'assets/*.xbm',
            'fonts/*.ttf',
        ],
    },
    entry_points={
    'console_scripts': [
        'chromacs=chromacs.chromacs_13c:main',
        'chromacs-overlap-expr=chromacs.chromacs_overlap_expr:main',
        ],
    },
)


