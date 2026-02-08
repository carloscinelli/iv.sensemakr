## Test environments

- local macOS (x86_64-apple-darwin20), R 4.4.1
- GitHub Actions: macOS-latest (R release), windows-latest (R release, R 4.1), ubuntu-latest (R devel, release, oldrel-1 through oldrel-4)

## R CMD check results

0 errors | 0 warnings | 4 notes

## Notes

- New submission.
- "unable to verify current time" reflects lack of internet access during local check.
- Markdown top-level file checks require pandoc; not available in this local environment.
- HTML validation issues are due to R's HTML generation, not package-specific problems.
