---
hide:
  - toc
---

# Flowchart to contribute to CABLE's documentation

```mermaid
   %% Create a graph for the legend of the main graph
   flowchart LR
      blank1[" "]
      uniq[action done <br> only once ever <br> at the start <br> to get the code locally]:::uniq;
      blank2[" "]

      blank1---uniq;
      subgraph Legend
         sevterm[action done <br> several times per issue <br> in a terminal/text editor]:::term;
         sevgit[action done <br> several times per issue <br> on GitHub]:::github;
         onceterm([action done <br> once per issue <br> in a terminal/text editor]):::term;
         oncegit([action done <br> once per issue <br> on GitHub]):::github;
         question[\question with multiple outcomes/]:::question;

         uniq --- sevterm --- sevgit --- question;
         uniq --- onceterm --- oncegit --- question;
      end
      question --- blank2

      classDef default fill:#FFFDE7, stroke:#FFF59D;
      classDef uniq fill:#D81B60, stroke:#880E4F, color:#FFFFFF;
      classDef github fill:#388E3C,stroke:#1B5E20, color:#FFFFFF;
      classDef term fill:#F44336, stroke:#B71C1C, color:#FFFFFF;
      classDef question fill:#6D4C41, stroke:#3E2723, color:#FFFFFF; 
;
```

```mermaid
   flowchart TB

      %% Define all the nodes first
      clone[Clone repository]:::uniq ;
      issue([Open an issue <br> to explain your work]):::github;
      branch([Create a branch]):::github;
      fetch([Fetch branch <br> to your local repository: <br> git fetch]):::term;
      locbranch(["Create local branch: <br> git checkout &lt;branchname&gt;"]):::term;
      edits[Make your edits]:::term;
      commit["Commit your edits: <br> git commit -a -m &lt;message&gt;"]:::term;
      push[Push your edits <br> git push]:::term;
      first[\Is it first push?/]:::question;
      pr([Create a pull request]):::github;
      preview[\Is preview correct?/]:::question;
      review[Request review]:::github;
      reviewchange[\Review asks for code changes <br> or has comments?/]:::question;
      revieweredits[\Reviewer has committed some edits? /]:::question;
      update[Pull updates to your local repository: <br> git pull]:::term;
      merge([Merge branch]):::github;
      updatemain([Pull update to your local repository: <br> git checkout main <br> git pull]):::term;

      %% Define the graphs using the nodes' names
      clone --> issue;
      merge --> updatemain;

      subgraph Review
         review --> reviewchange;
         reviewchange -->|No|merge;
      end

      subgraph Editing
         edits --> commit --> push --> first;
         first -->|Yes|pr --> preview;
         first -->|No|preview;
         preview -->|Yes|review;
         preview -->|No|edits;
         reviewchange -->|Yes|revieweredits;
         revieweredits -->|Yes|update --> edits;
         revieweredits -->|No|edits;
      end

      subgraph Setup
         issue --> branch --> fetch--> locbranch--> edits;
      end

      classDef default fill:#FFFDE7, stroke:#FFF59D;
      classDef uniq fill:#D81B60, stroke:#880E4F, color:#FFFFFF;
      classDef github fill:#388E3C,stroke:#1B5E20, color:#FFFFFF;
      classDef term fill:#F44336, stroke:#B71C1C, color:#FFFFFF;
      classDef question fill:#6D4C41, stroke:#3E2723, color:#FFFFFF; 

;
```

## Clone repository

![clone](../assets/clone.png)
Clone the repository to your local working space. Make sure to choose the appropriate protocol (HTTPS or SSH) for connecting to the remote repository. Note, you need to setup an access token to use HTTPS and SSH keys to use SSH.

## Open issue

![issue](../assets/issue.png)
Before starting new work, open an issue to explain what you are planning on working on. This avoid potential duplication of effort.

## Create a branch

![newbranch](../assets/newbranch.png)
From the issue screen on GitHub, you can directly create a branch for that issue from the option in the right-hand side menu. GitHub suggests a branch name you can edit, but please keep the issue number in the branch name.

## Fetch the branch locally

You have created a branch on the GitHub repository (aka the remote repository), you now need to update the repository on your work space with that information:

```bash
git fetch
```

## Create the local branch

You want the local branch to sync with the remote branch you have created previously. For this, use `git checkout` with the name of the branch you have created:

```bash
git checkout <branchname>
```

## Commit your edits

You need to record your edits in git, this is called `commit`. It is recommended to do this regularly as it gives some safety to reverse changes:

```bash
git commit -a -m "your commit message"
```

## Push your edits to the remote repository

It is recommended to push your changes back to the remote repository often as it provides a backup of the work, the ability to work from different computers and makes it easier to collaborate on some development:

```bash
git push
```

The first time you push back some work on a branch, consider opening a pull request. This allows potential collaborators or helpers to find your work easily.

![pr](../assets/pr.png)

You can update the pull request by simply pushing more commits to the same branch.

## Check the preview

The pull request will build a preview of your work *merged* with the main branch of the repository. Please check that your work is rendered correctly. Once a preview is ready, you will see the following comment in the pull request providing the path to the preview:

![preview](../assets/preview.png)

To preview changes to the API documentation (ie. documentation in the CABLE source code), you need to append `/api` to the path provided.

## Ask for review

![review](../assets/review.png)

Once you are satisfied with your work, ask for a review. By putting a submission in, you are responsible for being responsive to any comments or edit changes suggested by the reviewer. Remember your work will not be accepted into the main deployment branch until a reviewer approves it.

## Merge your branch

Once the reviewer(s) has(have) accepted your changes, you can merge your work into the main branch. Feel free to choose the merge method you prefer. "Merging" is the simplest and should be used if you don't understand what the other methods do.

## Update your local repository

Finally, don't forget to update your local repository to sync the main branch with the state of the remote repository. For this, you need to checkout main and then pull from the remote repository:

```bash
git checkout main
git pull
```
