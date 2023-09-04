from src.classical_ga import CodonOptimization
from src.codon_bqm import DWaveBQM, QiskitBQM, ibm_score
from src.constants import *
from Bio import SeqIO
import numpy as np

# In order to use D-Wave or IBM, you must have access to appropriate
# libraries/devices.
use_dwave = False
use_ibm = True
use_ga = True
use_paper_seqs = False

if use_paper_seqs:
    # Path to fasta
    fasta_path = 'sample_seqs/covid_sequences.fasta'
    seq_list = [str(bio_seq.seq) for bio_seq in SeqIO.parse(fasta_path, 'fasta')]
    seq_name = [str(bio_seq.id) for bio_seq in SeqIO.parse(fasta_path, 'fasta')]
else:
    seq_list = ['AGM', 'GSK']
    seq_name = ['test_1', 'test_2']

for p_seq, seq_name in zip(seq_list, seq_name):

    if use_dwave:
        # Must have D-Wave libraries installed. If using hybrid or QPU methods,
        # must run code on Leap interface with appropriate resource access.
        # Quantum code is programmed to run n_execs times
        dwave_results = DWaveBQM(p_seq, exact=True, hybrid=False)
        print('D-Wave:',dwave_results.score, dwave_results.score_mean, dwave_results.score_std)

    if use_ibm:
        # Must have qiskit installed or run code on IBM Experience.
        # QiskitBQM will run 1 time.
        ibm_results = QiskitBQM(p_seq, exact=False, noise=False)
        print('IBM:',ibm_results.score, ibm_results.exact_score)

    if use_ga:
        # Run classical GA n_execs times
        c_scores = [CodonOptimization(p_seq).score for _ in range(n_execs)]
        print('GA:',min(c_scores),np.mean(c_scores), np.std(c_scores))

