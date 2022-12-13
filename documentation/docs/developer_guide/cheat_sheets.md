# Cheat sheets

## Git

A minimum set of git commands used in the workflow for CABLE is given in the following table:

| Command | Description |
| ------- | ----------- |
| git clone <URL> | create a working copy of a remote repository |
| git checkout <branchname\> | Switch to a different branch <br> creates the branch if it does not exist |
| git add . | Add all modifications to the index before committing to git |
| git commit | Commit the index. You need to add an informative comment |
| git push | Send your local changes to the remote repository |
| git pull | Update your checkout branch with the state on the remote repository  |
| git branch -d <branchname\> | Delete a local branch. You need to have a different branch <br> currently checked out |

Additional commands and information can be found in [this "simple guide"][git-simple-guide]
## FORD

For details about the configuration or more explanation of the syntax, please refer to [the full documentation for FORD][ford-docs].

This cheat-sheet is intended for use for the documentation of the CABLE land surface model and focuses on the features needed for that purpose.

| Syntax | Description
| --- | --- |
| `!!` | **Comment line** for documentation <br> **after** the code to document |
| `!*` | **Block of comments** for documentation <br> **after** the code to document |
| `!>` | **Comment line** for documentation <br> **before** the code to document |
| `!|` | **Block of comments** for documentation <br> **before** the code to document |
| `\(...\)` | Inline mathematics notation |
| `\[...\]` | Mathematics notation on its own line |
| `\begin{equation}...\end{equation}` | Numbered equation on its own line |

## Markdown

For a quick reference to Markdown, please use [this cheat sheet][md-cheatsheet].

[md-cheatsheet]: https://www.markdownguide.org/cheat-sheet/
[ford-docs]: https://forddocs.readthedocs.io/en/latest/index.html
[git-simple-guide]: https://rogerdudler.github.io/git-guide/
