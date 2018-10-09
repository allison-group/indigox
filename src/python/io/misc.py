__all__ = ["LoadFile"]

## \brief Load a file.
#  \details Loads the given file line-by-line. Comments are removed from
#  the line.
#  \param path a Path object for the file to load.
#  \param comment character(s) to treat as comments. Anything after the
#  first appearance of these characters in a line will be ignored.
#  \param blanks include blank lines in the output. Set to False, this
#  will ignore lines with only whitespace, even if \p strip is False.
#  \param strip remove leading and trailing whitespace from the lines.
#  \return yields the loaded file, line by line.
def LoadFile(path, comment='#', blanks=False, strip=True):
  with path.open('r') as f:
    for line in f:
      if comment is not None and comment in line:
        line = line[:line.index(comment)]
      if strip:
        line = line.strip()
      if blanks:
        yield line
      elif line:
        yield line

