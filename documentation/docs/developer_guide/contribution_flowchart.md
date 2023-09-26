# Guidelines to contribute to CABLE

!!! warning "Pre-requisite"

    Before contributing to CABLE, please ensure you have followed all the steps to [setup Git and GitHub][git-training] given by ACCESS-NRI. Failing this, some of the commands described here may require additional steps or options.

These guidelines are written from the perspective of a new contributor to CABLE wanting to start their first development. Notes are added when guidelines for subsequent contributions differ.

:material-sim-alert: The CABLE documentation is an integral part of the CABLE code and its repository. The present guidelines apply for changes to the scientific code as well as to the documentation. All changes to the scientific code **must** be accompanied by adequate changes to the documentation. The changes to the documentation are expected to be part of the same set of changes as the scientific modifications. The documentation changes **must** not be submitted separately, except for corrections. You can refer to [the documentation guidelines][doc-guidelines] to know what level of documentation is expected and how to document code changes.

!!! info "Resources"

    Please refer to [this cheat sheet page][cheatsheet] for quick references to Git, Markdown and FORD syntax.

First, here is a flowchart explaining how the various steps of the workflow interact together. More details are provided for each step in the following sections of this page.

```mermaid
   flowchart TD

   Idea[Explain your work in an issue]
   Workspace[Create a place for your work for you in the repository]
   Work[Do your work, record it and check it]
   Review[Get a review on your work]
   Merge[Get your work into <br> the main development version of CABLE]
   FinalUpdate[Get ready for your next project]

   Idea --> Workspace --> Work;
   Work -->|Incorrect|Work;
   Work -->|Correct|Review;
   Review -->|Apply requested changes|Work;
   Review -->|Approved|Merge --> FinalUpdate;

   linkStyle 4 stroke:red,color:red;
   linkStyle 3 stroke:green,color:green;
   linkStyle 2 stroke:red,color:red;
   linkStyle 5 stroke:green,color:green;
;
```
!!! info "Legend of flowcharts"

    The legend for all subsequent flowcharts can be found in the [legend section](#legend-for-the-flowcharts)

## Explain your work

You first need to tell the community what you will be working on. For this you need to open an Issue on the CABLE's GitHub repository.

### Open an issue

<figure markdown>
  ![Image title](../assets/issue.png){ width="90%", align=right }
</figure>
Make sure to give a description that is detailed enough. The Issue should describe what you are trying to do and why and not how you are doing it. 

For example, if you want to implement an algorithm by Doe J. et al (1950), say so in your issue and ideally put a link to the publication (using a DOI). You can add why you think this is a good addition to CABLE: is it adding a completely new process? Is it adding a new parametrisation for an existing process?

You can cite some of the source code if it helps explain what the issue is.

!!! failure "Do not list the code changes you made to fix the issue"

    An issue should not contain the detail of the code changes made and the modified files. That information is contained in the Git history of the code and does not need to be replicated by hand.

    An issue should be used to discuss and plan the work before it is implemented so the work done to fix the issue should not be commented on here but in the pull request.

## A place for your work in the repository

```mermaid
   flowchart TB

      %% Define all the nodes first
      clone[Clone repository]:::uniq ;
      locbranch(["Create local branch: <br> git checkout &lt;branchname&gt;"]):::term;
      work[Ready to start work];

      clone --> locbranch--> work;

      classDef default fill:#FFFDE7, stroke:#FFF59D;
      classDef uniq fill:#D81B60,stroke:#880E4F,color:#FFFFFF;
      classDef github fill:#388E3C,stroke:#1B5E20,color:#FFFFFF;
      classDef term fill:#F44336,stroke:#B71C1C,color:#FFFFFF;
      classDef question fill:#6D4C41,stroke:#3E2723,color:#FFFFFF; 

      click work "https://cable.readthedocs.io/en/latest/developer_guide/contribution_flowchart/#do-and-record-your-work" "Doing work flowchart"
;
```

### Clone the repository

!!! note "First time only"

    You only need to clone the repository on your work machine the first time you work on a development for CABLE. The subsequent times, you can work from the same copy of the repository.

On [CABLE's GitHub main page][CABLE-repo], click the `code` green button, choose the SSH protocol and copy the URL you need:
<figure markdown>
  ![Image title](../assets/clone.png){ width="90%", align=right }
</figure>
!!! warning "Change the protocol after cloning"

    It is possible to change the access protocol to GitHub in your cloned repository if you realise you have chosen the wrong protocol, ie.:
    
    - you have cloned using the HTTPS protocol but have SSH keys setup with GitHub
    - or you have cloned using the SSH protocol but you have a personal access token setup with GitHub

    To do this, check our FAQs. 

On your local machine in a terminal, clone the repository:

```bash
git clone <URL provided>
```

### Create a branch for your work

You want to create your own branch for your work, that's your workspace within the shared repository. For this, in a terminal on your machine, use `git checkout` with the name of the branch you want:

```bash
git checkout <branchname>
```

???+ tip "Branch name"

    It is best to start the branch name with the issue number so the link between the two is obvious.

## Do, record and check your work

You are now ready to start implementing your modifications in the CABLE code and its documentation.

We recommend that you record your modifications on your machine with Git frequently as it provides with the ability to undo part of the changes easily. We also recommend to update your work on GitHub frequently. It ensures the community knows you are actively working on this feature, it provides you with a backup and cloud storage for your work, it enables you to check your modifications against some basic tests automatically during development.

```mermaid
   flowchart TB

      edits[Modify the source code]:::term;
      add["Add your edits to git's index: <br> git add ."]:::term;
      commit["Commit your edits locally: <br> git commit"]:::term;
      push[Push your edits to GitHub: <br> git push]:::term;
      pr([Create a pull request]):::github;
      tests[\Do the tests pass?/]:::question;
      review[Ready for review];

      edits --> add --> commit --> push;
      push --> |first time|pr --> tests;
      push ----> |subsequent|tests;
      tests -->|Yes|review;
      tests -->|No|edits;

      classDef default fill:#FFFDE7, stroke:#FFF59D;
      classDef uniq fill:#D81B60, stroke:#880E4F, color:#FFFFFF;
      classDef github fill:#388E3C,stroke:#1B5E20, color:#FFFFFF;
      classDef term fill:#F44336, stroke:#B71C1C, color:#FFFFFF;
      classDef question fill:#6D4C41, stroke:#3E2723, color:#FFFFFF; 

      linkStyle 7 stroke:red,color:red;
      linkStyle 6 stroke:green,color:darkgreen;

      click review "https://cable.readthedocs.io/en/latest/developer_guide/contribution_flowchart/#review" "Review process"

;
```

### Record your modifications

You need to record your edits in git, this is called `commit`. It is recommended to do this regularly as it gives some safety to reverse changes. The recommendation is to commit every logical part of your modifications. For example you could commit together: code to read in new input variables or files or the declaration of new variables.

In Git, committing is a 2-step process with a first step to add what you want to commit to the index and then commit what is in the index.

```bash
# To add all the modified files to the next commit
git add .
git commit
```

Then please add an informative message in the editor that opens automatically. You need to give some information on the nature of the changes. With the previous commit examples, it could be: "Read in variable X" or "Declare variables for process Y". You can give a short title to the commit followed by a longer description after an empty line.

!!! failure "Do not repeat what Git stores automatically"

    You do not need to list the following in your commit message: 

    - the date and time
    - the name of author
    - the email of the author
    - the files modified in the commit

!!! example "Commit message example"

    #1000- Add global variables for Doe J.'s equation

### Update the remote repository

It is recommended you update the remote repository with your changes often as it provides a backup of the work, the ability to work from different computers and makes it easier to collaborate on some development.

Updating the remote repository is called "push" for git and is simply done as is:

```bash
git push
```

### Open a pull request

!!! note "Once per branch only"

    You need to create the pull request only once per branch. Subsequent push from the same branch will update the status of the previously created pull request.

In addition to updating the remote repository with your changes, we recommend creating a pull request **at the start of your work**. In doing this, you provide yourself with:

- a place to ask for technical help on your development work
- a place to collaborate with others
- automated tests
- check if your changes conflict with the current development version of CABLE

Discussions about implementation details of a sub-feature should happen in a pull request rather than in the issue. This is because the pull request shows what has been done and comments can be left attached to specific lines of the code.

<figure markdown>
  ![Image title](../assets/pr.png){ width="90%", align=right }
  <figcaption>How to create a pull request</figcaption>
</figure>

### Check the preview

The pull request will build a preview of your work *merged* with the main branch of the repository. Please check that your work is rendered correctly. Once a preview is ready, you will see the following comment in the pull request providing the path to the preview:

<figure markdown>
  ![Image title](../assets/preview.png){ width="90%", align=right }
</figure>

To preview changes to the API documentation (ie. documentation in the CABLE source code), you need to append `/api` to the path provided.

## Review

All modifications to the documentation no matter how large or small need to be reviewed by another CABLE collaborator.

```mermaid
   flowchart TB

      %% Define all the nodes first
      review[Request review]:::github;
      update[Pull updates to your branch to your local repository: <br> git pull]:::term;
      merge([Merge branch]):::github;
      updatemain([Sync your local repository with the remote: <br> git checkout main <br> git pull]):::term;
      edits[Back to editing your branch];
      delete([Delete the branch on GitHub]):::github;

      %% Define the graphs using the nodes' names
      review --> |Need to put code changes|update --> edits;
      review -->|Approved|merge;
      merge --> delete --> updatemain;

      classDef default fill:#FFFDE7, stroke:#FFF59D;
      classDef uniq fill:#D81B60, stroke:#880E4F, color:#FFFFFF;
      classDef github fill:#388E3C,stroke:#1B5E20, color:#FFFFFF;
      classDef term fill:#F44336, stroke:#B71C1C, color:#FFFFFF;
      classDef question fill:#6D4C41, stroke:#3E2723, color:#FFFFFF; 

      linkStyle 0 stroke:red,color:red;
      linkStyle 2 stroke:green,color:darkgreen;
      linkStyle 1 stroke:red,color:red;
      linkStyle 3 stroke:green;
      linkStyle 3 stroke:green;

      click edits "https://cable.readthedocs.io/en/latest/developer_guide/contribution_flowchart/#do-and-record-your-work" "Doing work flowchart"

;
```

### Ask for review

<figure markdown>
  ![Image title](../assets/review.png){ width="90%", align=right }
</figure>

Once you are satisfied with your work, ask for a review. By putting a submission in, you are responsible for being responsive to any comments or edit changes suggested by the reviewer. Remember your work will not be accepted into the main deployment branch until a reviewer approves it.

The reviewer might make code changes suggestions. It is recommended to accept all the suggestions you agree with first and then, work on your local copy for any additional required modifications.

If you need to ask for a second review from the same reviewer, in the reviewer's list, click on the icon for this purpose.

### Merge your branch

Once the reviewer(s) has(have) accepted your changes, you can merge your work into the main branch. Feel free to choose the merge method you prefer. "Merging" is the simplest and should be used if you don't understand what the other methods do.

### Delete the branch

Once the branch has been merged with main, make sure to delete the branch and start any other piece of work from a different branch. Deleting the branch can be done from GitHub and the pull request screen.

### Update your local repository

Finally, don't forget to update your local repository to sync the main branch with the state of the remote repository. For this, you need to checkout main and then pull from the remote repository:

```bash
git checkout main
git pull
```

At this stage, you can also delete your local branch which has become unnecessary:

```bash
git branch --delete <branchname>
```

## Legend for the flowcharts

```mermaid
   %% Create a graph for the legend of the main graph
   flowchart LR
      uniq[action done <br> only once ever <br> at the start <br> to get the code locally]:::uniq;
      sevterm[action done <br> several times per issue <br> in a terminal/text editor]:::term;
      sevgit[action done <br> several times per issue <br> on GitHub]:::github;
      onceterm([action done <br> once per issue <br> in a terminal/text editor]):::term;
      oncegit([action done <br> once per issue <br> on GitHub]):::github;
      question[\question with multiple outcomes/]:::question;

      uniq --- sevterm --- sevgit --- question;
      uniq --- onceterm --- oncegit --- question;

      classDef default fill:#FFFDE7, stroke:#FFF59D;
      classDef uniq fill:#D81B60, stroke:#880E4F, color:#FFFFFF;
      classDef github fill:#388E3C,stroke:#1B5E20, color:#FFFFFF;
      classDef term fill:#F44336, stroke:#B71C1C, color:#FFFFFF;
      classDef question fill:#6D4C41, stroke:#3E2723, color:#FFFFFF; 
;
```


[git-training]: https://access-nri.github.io/Training/HowTos/GitAndGitHub/
[CABLE-repo]: https://github.com/CABLE-LSM/CABLE
[cheatsheet]: documentation_guidelines/cheat_sheets.md
[doc-guidelines]: documentation_guidelines/index.md
