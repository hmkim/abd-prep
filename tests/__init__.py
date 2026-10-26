"""Test suite for abd-prep.

Written with ``unittest`` so it runs with the standard library alone::

    python3 -m unittest discover -s tests -t . -v

``pytest`` collects the same classes if it is installed. No test needs the
network; the prototype regression test skips itself unless the original
RCSB files are present in ``~/.cache/abd-prep/``.
"""
