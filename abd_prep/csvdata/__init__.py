"""CSV (assay / results table) group -- M2.5. CLI group name: ``csv``.

The package is called ``csvdata`` because a directory named ``csv`` would
shadow the standard library module whenever the package directory lands
on ``sys.path`` (direct script execution, some test runners) -- the same
reason the PDB range parser is ``ranges`` and not ``select``.

* :mod:`schema`     -- standard columns, ``--map``, missing tokens, number classes (edit here for R8)
* :mod:`csvfile`    -- delimiter detection, header checks, row streaming (F19)
* :mod:`validate`   -- F20/F21 rules, ``clean`` writer, ``inspect``
* :mod:`crosscheck` -- F22 id/sequence cross-reference against a FASTA or PDB
"""
