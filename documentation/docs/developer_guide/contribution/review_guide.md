# Review guidelines

All contributions to the CABLE land surface model will be reviewed before inclusion in the model. The review process is intended to check that:

- [the coding standards][coding-standards] have been applied
- the documentation is understandable and follows at least the [minimum requirements][doc-min-req]
- [the required test results][tests-req] are provided
- the proposed changes address the problem explained in the issue and only this problem
- the proposed changes are correct
- the implementation of the changes follows the design of the CABLE model and will be maintainable

Once you have finished the implementation of your changes, please make sure:

- the description of the pull request is up-to-date. 
- all the required test results are either linked to in comments or copied in. If you have performed tests beyond the required tests, make sure these tests and their results are described in comments in your pull request.
- the automated checks pass. Contact the @CABLE-LSM/admin team if you need help understanding a failure of these tests.
- all conflicts have been solved. 

Once you are ready, [ask for a review by `CABLE-LSM/reviewers`][how-ask-review]. If you want a specific individual, you can use the same process to choose that person.

## Review process

### Review response

The reviewers team will try and reply to review requests quickly. If you think your review request has been lost, please ask for an update as a comment on your pull request and mention your reviewer (with @).

The review is likely to be an iterative process between the reviewer and the author. As such, it is appreciated if you, the author, are responsive once the review process starts. Make sure to keep discussions polite and courteous. You can reject suggestions from reviewers but the reasoning has to be explained and the rejection has to be agreed to by the reviewer.

[To incorporate code changes requested by the reviewer][how-review], you often need to incorporate these changes to your local repository and push them to GitHub again. It is possible for the reviewer to suggest changes that can be apply directly in GitHub. We recommend to:

- first apply the suggestions you agree with via GitHub
- then update your local branch with `git pull`
- apply other changes required by the review locally to your branch
- push the fully revised version to GitHub (`git push`)
- and finally [ask for a re-review][how-re-review] once you have resolved all points raised by the reviewer

To better understand the pull request interface for reviews on GitHub, please refer to [the GitHub documentation][github-review].

### Merge approved submission

Once a reviewer has approved the pull request, you can merge it. We prefer if the author merges the pull request as it provides you a last chance to spot an issue. However, the admin team will periodically merge any pull request approved some time ago and that has not been merged.

[coding-standards]: ../other_resources/coding_standards.md
[doc-min-req]: ../documentation_guidelines/index.md
[tests-req]: testing.md
[how-ask-review]: resources/how_to.md#ask-a-review
[how-checks]: resources/how_to.md#
[how-review]: resources/how_to.md#understand-a-review
[how-re-review]: resources/how_to.md#request-a-re-review
[github-review]: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests
