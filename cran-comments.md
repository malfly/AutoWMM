## Submission summary

This is the first submission of mypackage to CRAN.

## R CMD check results

There were no ERRORs, WARNINGs, or NOTEs.

## R-hub checks

I could not test using R-hub checks due to the following error:
- Error in curl::curl_fetch_memory(url, handle = handle) : 
  SSL peer certificate or SSH remote key was not OK: [builder.r-hub.io] SSL certificate problem: self signed certificate

The error persists when trying to use rhub::validate_email() directly.

## win-builder checks

No issues were found.

## CRAN comments

The examples for `countTree`, `wmmTree`, and `estTree` are computationally intensive and have been wrapped in `\donttest{}` or `\dontrun{}`. They should not be executed during CRAN checks.

## Special notes

None to report.
