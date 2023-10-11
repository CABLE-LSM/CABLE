# Plan your work

## Get your copy of the software

!!! note "First time only"

    You only need to clone the repository on your work machine the first time you work on a development for CABLE. The subsequent times, you can work from the same copy of the repository.

If you are a member of the CABLE-LSM organisation on GitHub, you can simply work from a locally cloned copy of the repository. You do not need to create a fork of the GitHub repository. To become a member of the CABLE-LSM organisation, please reply on [this issue][new_member]

??? tip "How to"
    [How to clone][how_to_clone]

## Open an issue to explain your work

Before starting any coding, you need to explain what you intend on doing in a [GitHub Issue][github_issues].

!!! warning "Check existing issues first"
    Before opening an issue, please check if a similar issue is opened. If you find an opened issue, do you need an additional issue or does your proposed work fits within the current existing issue? Can the existing issue be amended to include your idea while staying on topic and of a reasonable size? Feel free to discuss this in the existing issue if you are unsure.

The issue is where all discussions on planning how the work is going to be done should go. Comments in the issue can include discussions around:

* finding the best place(s) in the code to implement your proposed changes.
* identifying if the proposed work is of the appropriate size or if it should be split into smaller issues.
* for bug reports, identifying the best fix.
* if adding a new feature, finding the best name and location for the new namelist switch.

!!! failure "Do not list the code changes you made to fix the issue"

    An issue should not contain the detail of the code changes made and the modified files. That information is contained in the Git history of the code and does not need to be replicated by hand.

### Considerations

**Keep it on topic.** Issues should be self-contained logical units of work. Be wary of "also" which often indicates two unrelated changes are lumped together.

**Not too big and not too small.** Remember all units of work will be reviewed, try to keep these to reasonable sizes to avoid lengthy review delays.

**Explain the work.** If you want to propose a new algorithm, in your issue description give a reference(s) to the algorithm(s). Explain why you think if would be beneficial to have it in CABLE. Explain where you plan to implement it in the code if you know, or discuss this in the issue before starting working on it.

If the issue is about a bug or a code improvement, explain why it is a bug. Give the solution if you know it. Give a reproducible test case that highlights the issue if possible.

[how_to_clone]: resources/how_to.md#cloning-a-repository
[github_issues]: https://github.com/CABLE-LSM/CABLE/issues
[new_member]: https://github.com/CABLE-LSM/CABLE/issues/110
