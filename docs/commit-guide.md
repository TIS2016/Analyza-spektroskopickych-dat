### Commit guide
 
 
 [*"A commit message shows whether a developer is a good collaborator."* - Peter Hutterer](https://au.linkedin.com/in/peterhutterer/en)
 
In order to create a useful revision history, teams should first agree on a commit message convention.

---

**Basic guidlines are:**

1. Separate subject from body with a blank line
2. Limit the subject line to 50 characters
3. Capitalize the subject line
4. Do not end the subject line with a period
5. Use the imperative mood in the subject line
6. Wrap the body at 72 characters
7. Use the body to explain what and why vs. how

---

Number 1. is self-explanatory.
Number 2. is crucial (for seeing whole git message in one line).
Number 3. is something, I propose to alter.
Number 4. I agree with.
Number 5. explanation:
<dd>To remove any confusion, here's a simple rule to get it right every time.

A properly formed git commit subject line should always be able to complete the following sentence:

*If applied, this commit will* **your subject line here**

For example:
<ul><li>If applied, this commit will <b>refactor subsystem X for readability</b></li>
<li>If applied, this commit will <b>update getting started documentation</b></li>
<li>If applied, this commit will <b>remove deprecated methods</b></li>
<li>If applied, this commit will <b>release version 1.0.0</b></li></ul></dd>

Number 6. I agree with.
Number 7. I agree with.

---

Proposals of mine, are to add topic / file you are dealing with in your commit inside [ ] these brackets at the beggining of commit message. Include reference to issue you are dealing with at the end (with #no). And last one is to not capitalize begining of sentence after [ ] brackets.

#### Note:
**If you can not sumarize your commit with one simple sentence, maybe you should consider splitting your commit in multiple ones.**

Example of commit message:

```
[topic] summarize changes in 50 chars or  less

+ buffer utility class (eg. new implemented functionality, bug fix, etc.)
- removed feature X (eg. serializing data, etc.)

More detailed explanatory text, if necessary. Wrap it to about 72
characters or so.

See: #6
```

---

Explain the problem that the commit is solving. Focus on why you
are making this change as opposed to how (the code explains that).
Are there side effects or other unintuitive consequences of this
change?


Tell me your suggestions, or if you agree. Once we merge this into master, we will consider it as our go to source when commiting.

Sources used:
[Chris Beams - Git Commit](http://chris.beams.io/posts/git-commit/)


