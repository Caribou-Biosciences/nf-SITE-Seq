<!--
# Caribou-Biosciences/nf-SITE-Seq pull request

Many thanks for contributing to Caribou-Biosciences/nf-SITE-Seq!

Please fill in the appropriate checklist below (delete whatever is not relevant).
These are the most common things requested on pull requests (PRs).

Remember that PRs should be made against the dev branch, unless you're preparing a pipeline release.

Learn more about contributing: [CONTRIBUTING.md](https://github.com/Caribou-Biosciences/nf-SITE-Seq/tree/master/.github/CONTRIBUTING.md)
-->

## PR checklist

- [ ] This comment contains a description of changes (with reason).
- [ ] If you've added a new tool - have you followed the pipeline conventions in the [contribution docs](https://github.com/Caribou-Biosciences/nf-SITE-Seq/tree/main/.github/CONTRIBUTING.md)
- [ ] Ensure the test suite passes (`nextflow run . -profile test,docker --outdir <OUTDIR>`).
- [ ] Check for unexpected warnings in debug mode (`nextflow run . -profile debug,test,docker --outdir <OUTDIR>`).
- [ ] Usage Documentation in `docs/usage.md` is updated.
- [ ] Output Documentation in `docs/output.md` is updated.
- [ ] `CHANGELOG.md` is updated.
- [ ] `README.md` is updated (including new tool citations and authors/contributors).
