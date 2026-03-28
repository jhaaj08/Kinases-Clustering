Layer Probing Improves Kinase Functional Prediction with Protein Language Models
Ajit Kumar¹* and Indra Prakash Jha²*
¹Adobe, Sector 132, Noida 201304, Uttar Pradesh, India
²IIIT Delhi, New Delhi, India
*These authors contributed equally to this work.

Corresponding author: Ajit Kumar
Email: ajikumar@adobe.com
Tel: +918619190698

 
Author Biographical Notes
Ajit Kumar is a researcher at Adobe, Noida, India, working at the intersection of machine learning and computational biology, with a focus on protein representation learning and bioinformatics.
Indra Prakash Jha is a researcher at IIIT Delhi, India, working on protein function prediction and deep learning for biological sequence analysis.
 
Abstract
Accurate kinase functional classification from sequence alone remains a challenge, particularly for divergent or poorly characterized kinases. Protein language models (PLMs) such as ESM-2 offer rich sequence representations, but most downstream applications rely exclusively on final-layer embeddings — a choice rarely questioned despite evidence from NLP that intermediate layers capture distinct semantic information. Here, we systematically evaluate all 33 transformer layers of ESM-2 for kinase functional prediction, comparing unsupervised clustering and supervised classification across a homology-partitioned benchmark. Mid-to-late layers (20–33) outperform the final layer by 32% in unsupervised Adjusted Rand Index (ARI: 0.268 → 0.354) and improve homology-aware classification accuracy to 75.7%. We further integrate domain-aware sequence extraction and calibrated confidence estimation into a reproducible end-to-end pipeline. These findings challenge the default of final-layer reliance and provide a practical framework for PLM-based functional annotation, with implications for other protein families beyond kinases.

Keywords: Protein Language Models; Kinase Functional Classification; Layer Selection in Transformers; Domain-Specific Embeddings; ESM-2; Homology-Aware Evaluation
 
Key Points
•	Mid-to-late transformer layers (layers 20–33) of ESM-2 outperform the commonly used final layer by 32% in unsupervised clustering (ARI: 0.268 → 0.354) and improve supervised kinase functional classification accuracy to 75.7% under rigorous homology-aware evaluation.
•	Domain-level embedding extraction using Pfam-defined kinase domains substantially improves sequence representation quality compared to full-length protein embeddings, highlighting the importance of focusing on conserved catalytic cores.
•	Platt scaling reduces the Expected Calibration Error (ECE) from 0.154 to 0.110 (28% improvement), enabling reliable confidence-aware predictions where low-confidence sequences can be flagged for expert review.
•	Exemplar-based nearest-neighbour retrieval using mid-layer embeddings achieves a top-1 hit rate of 71.2% and MRR of 0.795, supporting interpretable and transparent functional annotation decisions.
•	The layer selection methodology presented here is model-agnostic and applicable to other transformer-based protein language models, providing a generalizable framework for improving PLM-based protein function prediction beyond kinases.
 
1. Introduction
1.1 Background
Protein kinases are critical enzymes for cell signaling and major drug targets, but predicting their specific function from the sequence remains challenging. Advances in protein language models (PLMs), particularly the Evolutionary Scale Modeling (ESM) family, have made it possible to extract biologically meaningful features directly from amino acid sequences. These models, trained on millions of sequences using masked language modeling objectives, have achieved state-of-the-art results in a range of protein prediction tasks. Yet, a common oversight is the over-reliance on the final transformer layer, which may not contain the most functionally relevant information.
1.2 The Layer Selection Problem
Most studies default to using the final transformer layer for downstream tasks, assuming it encapsulates the richest representation. However, transformer-based PLMs are hierarchical, and prior work in natural language processing (NLP) has shown that intermediate layers often encode more transferable or semantically relevant features [1–3]. This discrepancy raises a fundamental question: which transformer layers best capture biologically meaningful information for protein function prediction?
1.3 Objectives
Protein kinases play a central role in cellular signalling, and accurate functional classification is critical for understanding disease mechanisms and drug development. Although PLMs offer a powerful route to sequence-based functional inference, most existing applications rely solely on the final transformer layer - potentially overlooking biologically relevant information distributed across intermediate layers.

This study investigates the following questions:
1.	Can intermediate transformer layers in ESM-2 improve functional classification of kinase domains compared to the final layer?
2.	What is the optimal strategy for selecting and aggregating layer embeddings for kinase classification?
3.	How do different embedding strategies affect unsupervised clustering and supervised classification under homology-aware conditions?
To answer these, we build a reproducible pipeline incorporating domain-level extraction, layer probing, calibrated classification, and rigorous evaluation. Our goal is to provide a practical framework for improving functional predictions of protein kinases using PLM representations.
 
2. Related Work
2.1 Traditional Approaches to Kinase Classification
Previous methods have relied on homology-based annotation using tools like BLAST or HMMER, or motif-based heuristics curated from the literature. While these methods perform well for known families, they struggle with novel sequences and require significant expert intervention.
2.2 Deep Learning and PLM-Based Protein Function Prediction
With the rise of protein language models such as UniRep, TAPE, ProtBERT, and ProtT5, the focus has shifted to unsupervised representation learning. Embedding-based models have shown promise in capturing global and local sequence properties relevant to function. However, most applications simply extract the final-layer embedding - a design choice that may not be optimal.
2.3 Layer Selection in PLMs
Layer probing in natural language processing (e.g., BERT, GPT) has shown that semantic features often peak at intermediate layers. In protein ML, recent studies hint at similar trends, but systematic evaluations are rare. Our work contributes the first comprehensive analysis of layer-wise embedding utility in kinase classification.
 
3. Methods
3.1 Data Collection
We retrieved a curated dataset of protein kinase sequences from the UniProt SwissProt database (release October 2025) using the following query: reviewed:true AND (keyword:KW-0418 OR name:kinase*). Only canonical isoforms with sequence length greater than 100 amino acids were retained. Fragmentary sequences were excluded using UniProt flags. To focus on functionally relevant regions, Pfam domains PF00069 (Protein kinase domain) and PF07714 (Protein tyrosine kinase) were extracted using HMMER 3.3 with an E-value threshold of 0.001. CD-HIT 4.8.1 was used to remove redundancy at various identity thresholds (70%, 50%, 40%) to generate homology-aware splits. All tool versions and parameters were recorded for reproducibility.
3.2 Model Architecture
We used the ESM-2 650M model (esm2_t33_650M_UR50D), a 33-layer transformer encoder pretrained on UniRef50 using a masked language modeling (MLM) objective. Each residue is mapped to a 1280-dimensional vector. For sequences exceeding the model's maximum input length (1,022 residues), we applied a sliding window approach with overlap stitching. We compared different pooling strategies including mean over residues and [CLS] token extraction, and experimented with embeddings from different layer ranges: final layer (33), mid-layers (20–30), and mid-to-final layers (20–33).
Model Selection Rationale
We selected ESM-2 (650M parameters) [4] as our primary model for five reasons:
4.	State-of-the-art performance on protein tasks: ESM-2 achieves the highest accuracy among publicly available protein language models on CATH structure prediction (CATH 4.2: 87% top-1) and protein-protein interaction prediction. Meta AI's 2023 benchmark shows ESM-2 outperforms ESM-1b (+8%), ProtBERT (+12%), and ProtTrans (+6%) on functional annotation tasks.
5.	Evolutionary-scale training data: ESM-2 was trained on UniRef50 (2020) with ~50M sequences spanning diverse protein families, ensuring broad coverage of kinase evolutionary space.
6.	Appropriate architecture depth for layer probing: ESM-2's 33-layer transformer provides sufficient depth to explore intermediate representations, which is critical for our research question.
7.	Computational feasibility: The 650M parameter variant balances performance and accessibility. It runs on single consumer GPUs (e.g., NVIDIA RTX 3090/4090 with 24 GB VRAM) and processes ~20 sequences per minute.
8.	Established baseline for reproducibility: ESM-2 is the current de facto standard in protein ML (>2,000 citations in 2 years), with well-documented APIs and extensive community adoption.
Alternative models considered but not used: ESM-1b (33 layers, 650M, superseded by ESM-2); ProtBERT (12 layers, too shallow for layer probing); ProtTrans-XLNet-BFD (24 layers, less widely adopted); ESM-2 3B/15B (prohibitively expensive); AlphaFold2 embeddings (optimised for structure, not function); Ankh 2023 (promising but recent).
3.3 Training Details
For supervised classification, we trained a multinomial logistic regression model using scikit-learn with L2 regularization and balanced class weights. We adopted 5-fold stratified cross-validation on homology-aware train/test splits generated via CD-HIT clustering at 40% sequence identity. All splits and seeds were fixed for reproducibility. Calibration was performed using Platt scaling to adjust predicted probabilities, enabling reliability-aware deployment. Embeddings were standardized using StandardScaler before model training.
3.4 Evaluation Metrics
We used a comprehensive suite of evaluation metrics:
•	Unsupervised clustering: Adjusted Rand Index (ARI), Normalized Mutual Information (NMI), Purity, Hungarian Matching Accuracy, and Silhouette Score.
•	Supervised classification: Accuracy, Macro-F1, Weighted-F1, per-class precision/recall/F1, and top-3 accuracy.
•	Calibration: Expected Calibration Error (ECE) and log-loss before and after Platt scaling.
•	Exemplar retrieval: Mean Reciprocal Rank (MRR), top-k hit rate, and PR-AUC.
All metrics were computed using scikit-learn and statistically validated with bootstrapped confidence intervals and permutation tests as appropriate.
3.5 Clustering Setup
For unsupervised analysis, we applied standard K-Means clustering on the protein embeddings derived from different ESM-2 layer configurations. We fixed the number of clusters to 8, corresponding to the number of known kinase functional classes. All clustering was performed on length and domain-normalised embeddings, using cosine distance as the similarity metric. No ground-truth labels were used during clustering, and evaluation was performed post hoc using ARI, NMI, and Hungarian Matching Accuracy against ground-truth class labels.
 
4. Results
4.1 Intermediate Layers Improve Unsupervised Clustering
We evaluated unsupervised clustering performance using ESM-2 embeddings across different transformer layers. Averaging mid-to-late layers (layers 20–33) substantially improved clustering performance compared to using only the final layer. Specifically, the Adjusted Rand Index (ARI) increased from 0.268 (last layer only) to 0.354 (layers 20–33), a relative improvement of 32%. Normalized Mutual Information (NMI) and Hungarian Matching Accuracy also showed consistent gains.
•	Baseline (last layer only): Baseline (last layer only): ARI = 0.268, NMI = 0.360
•	Mid-to-late layers (20–33): Mid-to-late layers (20–33): ARI = 0.354, NMI = 0.501
Domain-level embeddings further improved separability over full-length sequences, highlighting the value of focusing on the conserved catalytic core.

 
Figure 1. Clustering performance (ARI) across ESM-2 layer selection strategies.
Alt text: Bar chart comparing Adjusted Rand Index (ARI) scores across four ESM-2 layer configurations: baseline (layer 33 only, ARI=0.268), all layers (ARI=0.312), mid-range layers 20-30 (ARI=0.353), and extended mid-to-late layers 20-33 (ARI=0.354). The extended mid-to-late configuration achieves the highest ARI, representing a 32% improvement over baseline.
4.2 Mid-Layer Averaging Boosts Supervised Classification
We trained a logistic regression classifier using different ESM-2 embedding strategies on homology-aware train/test splits (40% identity threshold). Using mid-to-late layer averages (20–33) yielded the highest accuracy and macro-F1:
•	Accuracy: 75.7%
•	Macro-F1: 0.668
•	Top-3 Accuracy: 94.8%
This configuration outperformed baselines such as k-NN (68.4%) and motif-only features (52.3%).

Table 1. Supervised classification performance across embedding strategies (40% sequence identity split). Bold values indicate best performance per metric.
Method	Accuracy	Macro-F1	Top-3 Accuracy
ESM-2 (Layer 33) + LR	70.2%	0.593	92.1%
ESM-2 (Layers 20–33) + LR	75.7%	0.668	94.8%
ESM-2 + k-NN (k=5)	68.4%	0.542	91.2%
Motif-only LR	52.3%	0.389	78.6%
LR = Logistic Regression; k-NN = k-Nearest Neighbours (k=5). All models evaluated under homology-aware train/test splits generated via CD-HIT at 40% sequence identity. Highlighted row (green) indicates the proposed method.

 
Figure 2. Confusion matrix for supervised classification across 8 kinase functional classes. Mid-layer averaged embeddings show high recall for most classes.
Alt text: Heatmap confusion matrix showing predicted versus true kinase functional class labels for 8 classes. Diagonal cells represent correct classifications. Most off-diagonal values are low, indicating high per-class recall. The matrix uses a colour gradient from white (zero) to dark blue (high count).
[Figure 3 — INSERT IMAGE HERE]
Figure 3. Classification performance across different homology identity thresholds (70%, 50%, 40%). The proposed approach generalises well across increasingly challenging splits.
Alt text: Line chart showing classification accuracy for ESM-2 mid-layer embeddings (layers 20-33) versus baseline (layer 33) across three homology-aware evaluation thresholds: 70%, 50%, and 40% sequence identity. Both metrics decrease at lower identity thresholds, but mid-layer embeddings consistently outperform the baseline at all thresholds.
4.3 Embedding Strategies and Pooling Comparison
We compared different pooling strategies and embedding sources:
•	Mean pooling: Mean pooling across residues performed best overall.
•	CLS token: CLS token was competitive for the final layer but underperformed for mid-layer embeddings.
•	Motif concatenation: Motif concatenation offered negligible gains (<2% ARI increase), indicating ESM-2 already captures these features.

[Figure 4 — INSERT IMAGE HERE]
Figure 4. Effect of pooling strategy on performance. Mean pooling consistently outperforms CLS token across both clustering and classification.
Alt text: Grouped bar chart comparing three pooling strategies (mean pooling, CLS token, motif concatenation) across two evaluation settings (clustering ARI and supervised accuracy). Mean pooling achieves the highest scores in both settings. CLS token performs similarly for the final layer but drops for mid-layer embeddings. Motif concatenation shows minimal improvement over the baseline.
4.4 Calibration Improves Decision Reliability
To ensure reliability in downstream applications, we applied Platt scaling to calibrate classification probabilities. This reduced the Expected Calibration Error (ECE) from 0.154 to 0.110 (28% improvement), and log-loss from 1.07 to 0.77. Approximately 18% of test sequences were flagged as low-confidence (probability < 0.7), enabling expert review.

[ 
Figure 5. Calibration curves before and after Platt scaling. Platt scaling reduces overconfidence and improves calibration.
Alt text: Reliability diagram showing two calibration curves: one before Platt scaling (dashed line, above the diagonal indicating overconfidence) and one after Platt scaling (solid line, closer to the perfect calibration diagonal). X-axis shows mean predicted probability and y-axis shows fraction of positives. Platt scaling reduces ECE from 0.154 to 0.110.
4.5 Exemplar Retrieval and Interpretability
We evaluated the embeddings using nearest-neighbour retrieval. Mid-layer embeddings (20–33) achieved:
•	Top-1 hit rate: Top-1 hit rate: 71.2%
•	Top-3 hit rate: Top-3 hit rate: 86.7%
•	MRR: MRR: 0.795

High cosine similarity (>0.992) reliably predicted family membership, suggesting that ESM-2 embeddings support interpretable and confident exemplar-based decisions.

[Figure 6 — INSERT IMAGE HERE]
Figure 6. Exemplar-based retrieval using cosine similarity for mid-layer embeddings (layers 20–33) shows high interpretability and precision.
Alt text: Precision-recall curve for exemplar-based nearest-neighbour retrieval using ESM-2 mid-layer embeddings (layers 20-33). The curve shows high area under the curve (PR-AUC), with top-1 hit rate of 71.2% and top-3 hit rate of 86.7%. A threshold line at cosine similarity 0.992 is shown, above which predictions are considered high-confidence.
 
5. Mathematical Formalization of Layer Averaging Strategy
5.1 Layer Selection and Averaging Framework
Let L ⊆ {1, 2, …, 33} be the set of transformer layers selected for averaging. We define four configurations:
•	L_baseline = {33} (last layer only)
•	L_all = {1, 2, …, 33} (all layers)
•	L_mid-range = {20, 21, …, 30} (mid-to-late layers)
•	L_extended = {20, 21, …, 33} (extended mid-to-late layers)

We seek to identify the optimal configuration that maximizes clustering performance: L* = argmax_L ARI(L).
5.2 Multi-Layer Embedding Extraction
Given a protein sequence x = (x₁, x₂, …, x_L) of length L, ESM-2 provides per-residue embeddings H⁽ℓ⁾ ∈ R^(L×d) at each layer ℓ ∈ L, where d = 1280. The layer-averaged embedding is computed as:

H̅ = (1/|L|) Σ_{ℓ∈L} H⁽ℓ⁾

For each residue i ∈ {1, …, L}, we define the averaged representation: h̅_i = (1/|L|) Σ_{ℓ∈L} h_i^(ℓ).
5.3 Sequence-Level Pooling
To obtain the final sequence embedding z ∈ R^d, mean pooling is applied across residues: z = (1/L) Σ_{i=1}^{L} h̅_i. This operation yields a length-invariant sequence representation that captures global contextual information.
5.4 Handling Long Sequences
For sequences exceeding the ESM-2 token limit (L_max = 1022), the input sequence is divided into W overlapping windows with stride s = 900 residues. The final sequence embedding is computed via length-weighted averaging: z_final = (Σ_{w=1}^{W} n_w · z^(w)) / (Σ_{w=1}^{W} n_w). Weighting by window length ensures that longer segments contribute proportionally to the final embedding, avoiding bias toward overlapping regions.
5.5 Variance Reduction and Statistical Justification
Averaging across multiple layers reduces embedding variance. Assuming independence across layer representations: Var(h̅_i) = (1/k²) Σ_{ℓ∈L} Var(h_i^(ℓ)) = σ²/k, where k = |L| and σ² is the variance of a single-layer embedding. Thus, layer averaging reduces variance by a factor of k. By the Central Limit Theorem (CLT), averaged embeddings converge to a Gaussian distribution with reduced variance, improving stability and robustness.
5.6 Empirical Results Summary
The following summarises the empirical clustering performance (ARI) for each layer configuration:
•	ARI(L_baseline) = 0.268
•	ARI(L_all) = 0.312
•	ARI(L_mid-range) = 0.353
•	ARI(L_extended) = 0.354

Conclusion: The optimal configuration is L* = L_extended = {20, 21, …, 33}, achieving a 32% improvement over the baseline (p < 0.001, permutation test, Cohen's d = 1.87).
[Figure 7 — INSERT IMAGE HERE]
Figure 7. Layer-wise ARI heatmap showing performance across all 33 ESM-2 layers and layer averaging configurations.
Alt text: Heatmap showing ARI scores for individual ESM-2 layers (1-33) on the x-axis and different layer averaging window configurations on the y-axis. Warmer colours indicate higher ARI. The mid-to-late layer range (20-33) shows consistently the highest ARI values, with a clear gradient from lower ARI in early layers to higher ARI in mid-to-late layers.
 
6. Conclusion
In this study, we systematically explored the influence of transformer layer selection on the performance of protein language model (PLM) embeddings for kinase functional classification. Our findings reveal that the commonly used final layer is not the most informative, and instead, embeddings derived from mid-to-late layers (specifically layers 20–33 of ESM-2) significantly enhance both unsupervised and supervised tasks.

By combining layer-wise averaging with calibrated classification and domain-aware embedding extraction, we present a practical and reproducible pipeline that outperforms traditional motif-based and single-layer approaches. The 32% gain in ARI and 6% gain in Macro-F1 over the final-layer baseline establish the value of probing PLM depth.

Importantly, our approach remains accessible, as it uses the 650M parameter ESM-2 model — suitable for execution on single-GPU systems — and generalisable, with methods applicable to other PLMs like ProtT5, ProtTrans, and Ankh.

Beyond performance improvements, our work highlights the need for reliable calibration and interpretability in protein ML workflows. The use of length-weighted sliding windows, confidence estimation, and exemplar-based retrieval strengthens the biological relevance of our predictions.

Overall, this study contributes methodological insights and practical tools to the protein function prediction community, providing a foundation for future work on probing, fine-tuning, and transferring knowledge across deep protein models.
 
Acknowledgements
We thank the contributors of UniProt, Pfam, and Meta AI's ESM repository for providing high-quality data and tools that enabled this research. Computational resources were supported by local GPU clusters at our institution.
 
Declarations
Funding: Not applicable.

Conflict of interest: The authors declare no competing interests.

Ethics approval: Not applicable.

Consent to participate: Not applicable.

Consent for publication: All authors consent to the publication of this work.

AI/LLM use disclosure: [INSERT: Disclose any use of AI tools in writing, coding, or analysis, or state "No AI tools were used in the preparation of this manuscript."]
 
Data and Code Availability
Data availability: Benchmark datasets, preprocessed kinase domain files, and all evaluation splits are available at: https://github.com/jhaaj08/Kinases-Clustering. Archived copy: Zenodo DOI: 10.5281/zenodo.17370925.

Code availability: Full training and evaluation pipeline (including Snakemake scripts, environment files, and figure generation) is available at: https://github.com/jhaaj08/Kinases-Clustering. Archived copy: Zenodo DOI: 10.5281/zenodo.17370925.
 
[1] Peters, M. E., Neumann, M., Iyyer, M., Gardner, M., Clark, C., Lee, K., & Zettlemoyer, L. (2018). Deep contextualized word representations. Proceedings of NAACL-HLT 2018, 2227–2237.
[2] Jawahar, G., Sagot, B., & Seddah, D. (2019). What does BERT learn about the structure of language? ACL Workshop on Deep Learning and Formal Languages.
[3] Rogers, A., Kovaleva, O., & Rumshisky, A. (2020). A primer in BERTology: What we know about how BERT works. Transactions of the Association for Computational Linguistics (TACL), 8, 842–866.
[4] Lin, Z., Akin, H., Rao, R., Hie, B., Zhu, Z., Lu, W., et al. (2023). Evolutionary-scale prediction of atomic-level protein structure with a language model. Science, 379(6637), 1123–1130.
[5] Meta AI. (2023). ESM: Evolutionary Scale Modeling. Available at: https://github.com/facebookresearch/esm
[6] The UniProt Consortium. (2023). UniProt: the Universal Protein Knowledgebase in 2023. Nucleic Acids Research, 51(D1), D523–D531.
[7] Mistry, J., Chuguransky, S., Williams, L., et al. (2021). Pfam: The protein families database in 2021. Nucleic Acids Research, 49(D1), D412–D419.
[8] Fu, L., Niu, B., Zhu, Z., Wu, S., & Li, W. (2012). CD-HIT: accelerated for clustering the next-generation sequencing data. Bioinformatics, 28(23), 3150–3152.
[9] Pedregosa, F., Varoquaux, G., Gramfort, A., et al. (2011). Scikit-learn: Machine Learning in Python. Journal of Machine Learning Research, 12, 2825–2830.
