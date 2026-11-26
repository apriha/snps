# Contributing

Contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit
will always be given.

## Bug reports

When [reporting a bug](https://github.com/apriha/snps/issues) please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

## Documentation improvements

`snps` could always use more documentation, whether as part of the official `snps` docs, in
docstrings, or even on the web in blog posts, articles, and such. See below for info on how to
generate documentation.

## Feature requests and feedback

The best way to send feedback is to file an issue at https://github.com/apriha/snps/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that code contributions are welcome :)

## Development

To set up `snps` for local development:

1. Fork [snps](https://github.com/apriha/snps) (look for the "Fork" button).

2. Clone your fork locally:

   ```bash
   $ git clone git@github.com:your_name_here/snps.git
   ```

3. Create a branch for local development from the `main` branch:

   ```bash
   $ cd snps
   $ git checkout main
   $ git checkout -b name-of-your-bugfix-or-feature
   ```

4. Install [uv](https://docs.astral.sh/uv/) (see the
   [installation guide](https://docs.astral.sh/uv/getting-started/installation/)) if you don't
   already have it, then set up a development environment and install the
   [pre-commit](https://pre-commit.com) hooks:

   ```bash
   $ uv sync
   $ uv run pre-commit install
   ```

   `uv sync` creates a virtual environment and installs `snps` (editable) plus the development
   dependencies. `pre-commit install` registers the git hook so linting and formatting run
   automatically on every commit.

5. When you're done making changes, run all the tests with:

   ```bash
   $ uv run pytest --cov-report=html --cov=snps tests README.md
   ```

   > **Note:** Downloads during tests are disabled by default. To enable downloads, set the
   > environment variable `DOWNLOADS_ENABLED=true`.

   > **Note:** If you receive errors when running the tests, you may need to specify the temporary
   > directory with an environment variable, e.g., `TMPDIR="/path/to/tmp/dir"`.

   > **Note:** After running the tests, a coverage report can be viewed by opening
   > `htmlcov/index.html` in a browser.

6. Lint, format, and check types. The pre-commit hook runs the linting and formatting automatically
   on commit; you can also run them on demand:

   ```bash
   $ uv run pre-commit run --all-files   # ruff lint, formatting, and file checks
   $ uv run mypy                         # type checking
   ```

7. Commit your changes and push your branch to GitHub:

   ```bash
   $ git add .
   $ git commit -m "Your detailed description of your changes."
   $ git push origin name-of-your-bugfix-or-feature
   ```

8. Submit a pull request through the GitHub website.

### Pull request guidelines

If you need some code review or feedback while you're developing the code, just make the pull
request.

For merging, you should:

1. Ensure tests pass.
2. Update documentation when there's new API, functionality, etc.
3. Add yourself to `CONTRIBUTORS.md` if you'd like.

## Adding or changing a reader

Each file format has a module in `src/snps/io/readers/`. `livingdna.py` is a complete example;
`ancestry.py` also adjusts the table after reading it.

To add a format:

1. Create `src/snps/io/readers/<vendor>.py` with a module docstring naming the vendor and its URL,
   and define `FORMAT = ReaderFormat(name=..., detect=..., parse=...)`:
   - `name`: the format's name, reported as the file's source unless the optional `source(header)`
     returns another.
   - `detect(header)`: returns `True` if `header.first_line` or `header.comments` (the file's
     leading comment lines) mark the file as this format.
   - `parse`: if the data is one delimited table, describe it with a `TableLayout` (its fields are
     documented in `_table.py`) and use `partial(parse_table, LAYOUT)`; if the table needs
     adjusting, write a function that calls `read_table(parse_input.source, LAYOUT)` and adjusts
     the result. A format that is not a table provides its own parser (see `vcf.py`).
2. Add `FORMAT` to `_FORMATS` in `src/snps/io/readers/__init__.py`. The first format whose `detect`
   matches wins, so put it before any format that would also match its files.
3. Add a small invented file to `tests/input/` and a `ReadCase` for it to `CASES` in
   `tests/support/cases.py`.
4. Add the name to `FORMAT_NAMES` in `tests/io/test_readers.py`. If another format's `detect` also
   matches the new header, add a case to `test_reader_precedence` there.
5. Add the vendor to "Supported Genotype Files" in `README.md`.

To change a reader, edit its layout or its `parse`, and add a `ReadCase` (with a new file in
`tests/input/` if needed) that shows the new behavior.

## Documentation

After the development environment has been setup, documentation can be generated via the
following command (the docs dependencies are listed in `docs/requirements.txt`):

```bash
$ uv run --with-requirements docs/requirements.txt sphinx-build -T -E -D language=en docs docs/_build
```

Then, the documentation can be viewed by opening `docs/_build/index.html` in a browser.
