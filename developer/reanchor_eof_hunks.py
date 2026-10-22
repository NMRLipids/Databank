"""
Re-anchor an end-of-file insertion so a review suggestion can carry it.

``autocomplete_expsim_metadata.py`` appends its generated ``bioschema_properties``
block to the end of a ``README.yaml``. ``git diff`` expresses an append as a pure
insertion -- ``@@ -41,0 +42,74 @@`` -- and reviewdog anchors an inserted hunk at
the line *after* the last unchanged one, which for an append is one line past the
end of the file. No review comment can attach to a line that is not in the pull
request's diff, so ``reviewdog/action-suggester`` dropped every suggestion
carrying the block while the job still reported success.

This filter rewrites such a hunk into a replacement of the file's last line: the
trailing context line becomes a deletion followed by an identical addition, so
reviewdog anchors on a line that exists and the suggestion reads as that line
followed by the block. Applying it gives back, byte for byte, the file the script
wrote -- only the way the change is expressed differs, never the result.

Hunk headers are left alone: turning one context line into a ``-``/``+`` pair
spends one old line and one new line exactly as the context line did, so both
counts in ``@@ -a,b +c,d @@`` stay correct.

Only the last hunk of each file is considered, because only there can added lines
reach the end of the file: git emits trailing context after every other change.
A hunk whose trailing additions already follow a ``-`` line is left as it is --
reviewdog anchors that on the deleted line, which is a real one.

Reads a unified diff on stdin and writes one on stdout::

    git diff | python developer/reanchor_eof_hunks.py | reviewdog -f=diff ...
"""

import sys


def reanchor_hunk(body):
    """One hunk's lines, with a trailing insertion re-anchored onto the line above.

    Returned unchanged when there is nothing to do: no trailing run of ``+``
    lines, no line above them to anchor on, or a run that already continues a
    ``-``/``+`` change group. A ``\\ No newline at end of file`` marker sits
    between the additions and whatever precedes them and so stops the scan,
    which is the wanted behaviour -- both spellings of a missing final newline
    already produce an anchor on a real line.
    """
    start = len(body)
    while start > 0 and body[start - 1].startswith("+"):
        start -= 1
    if start in (0, len(body)):
        return body
    anchor = body[start - 1]
    if not anchor.startswith(" "):
        return body
    text = anchor[1:]
    return body[:start - 1] + ["-" + text, "+" + text] + body[start:]


def reanchor(diff_text):
    """The diff with every file's end-of-file insertion re-anchored."""
    out, hunk_start = [], None

    def close(is_last_hunk_of_file):
        nonlocal hunk_start
        if hunk_start is not None and is_last_hunk_of_file:
            out[hunk_start:] = reanchor_hunk(out[hunk_start:])
        hunk_start = None

    for line in diff_text.splitlines(keepends=True):
        if line.startswith("diff --git "):
            close(True)
        elif line.startswith("@@"):
            # Another hunk follows in this file, so the one just read cannot be
            # the one that reaches the end of it.
            close(False)
            out.append(line)
            hunk_start = len(out)
            continue
        out.append(line)
    close(True)
    return "".join(out)


def main():
    sys.stdout.write(reanchor(sys.stdin.read()))


if __name__ == "__main__":
    main()
