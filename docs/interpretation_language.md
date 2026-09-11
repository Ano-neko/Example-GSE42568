# Interpretation language

## Allowed

- “Higher/lower expression in breast cancer than in normal breast tissue in this cohort.”
- “Adjusted cancer-minus-normal association after including the processing-run proxy.”
- “The gene set is enriched toward the cancer-up/normal-up end of the P1 moderated-t ranking.”
- “The result is stable across the named sensitivity models under the frozen operational criteria.”
- “The pattern is consistent with, but does not prove, a biological program described in the cited literature.”

## Not allowed without new evidence

- “Healthy controls,” “matched normal,” or “adjacent normal.”
- “Cancer causes this gene/pathway change.”
- “The pathway is activated/inhibited.”
- “This gene diagnoses breast cancer, predicts survival, or predicts treatment response.”
- “This result identifies a drug target” or any treatment recommendation.
- “Batch effects were removed.” Run was modelled; balance was not randomized and the inferential matrix was not ComBat-adjusted.
- Treating leading-edge genes as an independently FDR-controlled biomarker list.

All public prose must name the contrast direction and preserve the bulk-tissue, normal-provenance, run-confounding and processed-only QC limitations.
