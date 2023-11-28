# Review guidelines

All contributions to the CABLE land surface model will be reviewed before inclusion in the model. The review process is intended to check:

- [the coding standards][coding-standards] have been applied
- the documentation is understandable and follows at least the [minimum requirements][doc-min-req]
- [the required test results][tests-req] are provided
- the proposed changes address the problem explained in the issue and only this problem
- the proposed changes are correct
- the implementation of the changes follows the design of the CABLE model and will be maintainable

Once you have finished the implementation of your changes, please make sure the description of the pull request is up-to-date. Ensure all the required test results are either linked to in comments or copied in. If you have performed tests beyond the required tests, make sure these tests and their results are described in comments in your pull request.

Once you are ready, [ask for a review by `CABLE-LSM/reviewers`][how-review]. If you want a specific individual who isn't part of the reviewers team to review your pull request, you can use the same process. However, a following review from the reviewers team will be necessary to accept the changes.

## Review process

The reviewers team will try and reply to review requests quickly. If you think your review request has been lost, please ask for an update as a comment on your pull request and mention your reviewer (with @).

The review is likely to be an iterative process between the reviewer and the reviewee. As such, it would be appreciated if the reviewee is responsive once the review process starts. Make sure to keep discussions polite and courteous. Reviewees can reject suggestions from reviewers but the reasoning has to be explained and the rejection has to be agreed to by the reviewer.

Once a reviewer has given you a written agreement the pull request is ready to be merged, you can merge it. We prefer if the author of the code changes merges the pull request as that leaves them a last chance to change their mind if they discover a reason not to go ahead with the changes. However, the admins team will periodically merge any approved pull request that has not been merged.

[coding-standards]: ../other_resources/coding_standards.md
[doc-min-req]: ../documentation_guidelines/index.md
[tests-req]: testing.md
[how-review]: resources/how_to.md#ask-a-review
