# Contributing to __pe*ga*susQC__

We welcome contributions to __pe*ga*susQC__, in the form of issues, bug fixes, documentation or suggestions for enhancements. This document sets out our guidelines and best practices for such contributions.

It's based on the [Contributing to Open Source Projects Guide](https://contribution-guide-org.readthedocs.io/).

__pe*ga*susQC__ has the following modes of contribution:

- GitHub Commit Access
- GitHub Pull Requests

## Code of Conduct

Contributors to this project are expected to act respectfully toward others in accordance with the [Code of Conduct](https://github.com/GeoscienceAustralia/pegasusQC/blob/main/CODE_OF_CONDUCT.md).

## Submitting Bugs

### Due Diligence

Before submitting a bug, please do the following:

- Perform basic troubleshooting steps:
	- Make sure you're on the latest version. If you're not on the most recent version, your problem may have been solved already! Upgrading is always the best first step.
	- [__Search the issue
      tracker__](https://github.com/GeoscienceAustralia/pegasusQC/issues) to make sure it's not a known issue.

### What to put in your bug report

Make sure your report gets the attention it deserves: bug reports with missing information may be ignored or punted back to you, delaying a fix. The below constitutes a bare minimum; more info is almost always better:

- What version of Python are you using? For example, are you using Python 2.7, Python 3.7, PyPy 2.0?
- What operating system are you using? Windows (7, 8, 10, 32-bit, 64-bit), Mac OS X, (10.7.4, 10.9.0), GNU/Linux (which distribution, which version?) Again, more detail is better.
- Which version or versions of the software are you using? Ideally, you've followed the advice above and are on the latest version, but please confirm this.
- How can the we recreate your problem? Imagine that we have never used __pe*ga*susQC__ before and have downloaded it for the first time. Exactly what steps do we need to take to reproduce your problem?

## Contributions and Licensing

### Contributor License Agreement

Your contribution will be under our [license](https://github.com/GeoscienceAustralia/pegasusQC?tab=MIT-1-ov-file#readme) as per [GitHub's terms of service](https://docs.github.com/en/site-policy/github-terms/github-terms-of-service#6-contributions-under-repository-license).

### GitHub Commit Access

Commit access is limited to the project team.

### GitHub Pull Requests

- Pull requests may include copyright in the source code header by the contributor if the contribution is significant or the contributor wants to claim copyright on their contribution.
- All contributors shall be listed on the __pe*ga*susQC__ github repository.
- Unclaimed copyright, by default, is assigned to the main copyright holders as specified in the [license](https://github.com/GeoscienceAustralia/pegasusQC/blob/main/LICENSE).

### Version Control Branching

- Always make a new branch for your work, no matter how small. This makes it easy for others to take just that one set of changes from your repository, in case you have multiple unrelated changes floating around.
	- Don't submit unrelated changes in the same branch/pull request! If it is not possible to review your changes quickly and easily, we may reject your request.
- Base your new branch off of the appropriate branch on the main repository:
	- In general the released version of __pe*ga*susQC__ is based on the master (default) branch whereas development work is done under other non-default branches. Unless you are sure that your issue affects a non-default branch, base your branch off the master one.
- Note that depending on how long it takes for the dev team to merge your patch, the copy of master you worked off of may get out of date!
	- If you find yourself 'bumping' a pull request that's been sidelined for a while, make sure you rebase or merge to latest master to ensure a speedier resolution.

### Documentation

- documentation is managed in docs/, in markdown format, with tutorials in Jupyter-lab notebooks.
- Sphinx is used to generate the documentation.

### Code Formatting

- Please follow the coding conventions and style used in the __pe*ga*susQC__ repository.
- __pe*ga*susQC__ follows the PEP-8 guidelines
- 80 characters
- spaces, not tabs
- __pe*ga*susQC__, or pegasusQC instead of PegasusQC, pegasusqc, etc.

## Suggesting Enhancements

We welcome suggestions for enhancements, but reserve the right to reject them if they do not follow future plans for __pe*ga*susQC__.