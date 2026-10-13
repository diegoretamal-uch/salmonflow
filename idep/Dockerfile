# ──────────────────────────────────────────────────────────────
# SalmonFlow's iDEP: gexijin/idep 2.4.2 with a patched idepGolem
#
# The patch (https://github.com/diegoretamal-uch/idepGolem, branch
# salmonflow-full-deseq2) adds a "Full DESeq2 results" download to the
# DEG tab: every DESeq2::results() column per comparison, including
# lfcSE and the Wald stat used to rank genes for GSEA. Everything else
# is upstream iDEP 2.4.2, unchanged.
#
# Only the idepGolem package is replaced. The base image (R, DESeq2 and
# the bundled annotation database) is reused as-is, so building this
# takes about a minute once the base is present:
#   docker pull gexijin/idep:latest     # one time, large
#   docker compose build idep           # or: docker build -t salmonflow-idep idep/
#
# The base is pinned by digest because the fork is branched from the
# idepGolem commit this exact image was built from (v240_deployed,
# f442289). A newer base would need the fork rebased to match.
# ──────────────────────────────────────────────────────────────
FROM gexijin/idep@sha256:ea5b1df64ed7ab1f12057132c08dc934bc0ee74e0b863a4ecce4517dc5a0d847

ARG IDEPGOLEM_REPO=diegoretamal-uch/idepGolem
ARG IDEPGOLEM_SHA=49fc63e30275340592c18950e83bba3110ddfbad

# dependencies = FALSE: everything idepGolem needs is already in the
# base image; letting remotes resolve them would upgrade packages and
# change iDEP's numbers.
RUN R -q -e "remotes::install_github('${IDEPGOLEM_REPO}@${IDEPGOLEM_SHA}', \
      dependencies = FALSE, upgrade = 'never', force = TRUE)" \
 && R -q -e "stopifnot(packageDescription('idepGolem')\$RemoteSha == '${IDEPGOLEM_SHA}')"
