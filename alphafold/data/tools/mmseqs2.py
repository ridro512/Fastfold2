"""MMseqs2 wrapper to replace HHblits in AlphaFold MSA generation."""

import glob
import os
import subprocess
import tempfile
from typing import Sequence, Optional, Mapping, Any
from alphafold.data.tools import utils


class Mmseqs2HhblitsReplacement:
    """MMseqs2 runner specifically designed to replace HHblits searches."""
    
    def __init__(self, 
                 binary_path: str,
                 databases: Sequence[str],
                 n_cpu: int = 8,
                 sensitivity: float = 8.0,
                 max_sequences: int = 10000,
                 use_gpu: bool = False):
        """Initialize MMseqs2 runner."""
        self.binary_path = binary_path
        self.databases = databases
        self.n_cpu = n_cpu
        self.sensitivity = sensitivity
        self.max_sequences = max_sequences
        self.use_gpu = use_gpu
        
        # Verify binary exists
        if not os.path.exists(self.binary_path):
            raise FileNotFoundError(f"MMseqs2 binary not found at {self.binary_path}")
        
        # Verify databases exist
        valid_databases = []
        for db in self.databases:
            if os.path.exists(db):
                valid_databases.append(db)
            else:
                print(f"Warning: Database not found at {db}")
        
        if not valid_databases:
            raise ValueError("No valid databases found")
        
        self.databases = valid_databases
    
    def query(self, input_fasta_path: str) -> Mapping[str, Any]:
        """Run MMseqs2 search against databases.
        
        Args:
            input_fasta_path: Path to input FASTA file
            
        Returns:
            Dictionary with 'a3m' key containing the A3M format MSA
        """
        # Create temporary directory for MMseqs2 operations
        with tempfile.TemporaryDirectory() as tmp_dir:
            all_results = []
            
            for i, database in enumerate(self.databases):
                print(f"Searching database {i+1}/{len(self.databases)}: {os.path.basename(database)}")
                
                # Create temporary files for this database search
                result_prefix = os.path.join(tmp_dir, f'result_{i}')
                result_m8 = f"{result_prefix}.m8"
                tmp_search_dir = os.path.join(tmp_dir, f'tmp_{i}')
                
                # MMseqs2 easy-search command
                cmd = [
                    self.binary_path, 'easy-search',
                    input_fasta_path,
                    database,
                    result_m8,
                    tmp_search_dir,
                    '--threads', str(self.n_cpu),
                    '-s', str(self.sensitivity),
                    '--max-seqs', str(self.max_sequences),
                    # FIXED: Use qaln,taln for aligned sequences
                    '--format-output', 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln',
                    '-e', '1e-3',  # E-value threshold
                    '--comp-bias-corr', '1',
                    '--remove-tmp-files',  # Clean up temporary files
                    '--alignment-mode', '3',  # Better alignment quality
                    '--min-seq-id', '0.0'     # Don't filter by sequence identity
                ]
                
                # Add GPU flag if requested
                if self.use_gpu:
                    cmd.extend(['--gpu', '1'])
                
                try:
                    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
                    if os.path.exists(result_m8) and os.path.getsize(result_m8) > 0:
                        all_results.append(result_m8)
                        print(f"  Found {self._count_hits(result_m8)} hits")
                    else:
                        print(f"  No hits found")
                except subprocess.CalledProcessError as e:
                    print(f"Warning: MMseqs2 search failed for database {database}")
                    print(f"Error: {e.stderr}")
                    continue
            
            # Convert results to A3M format
            a3m_content = self._convert_to_a3m(all_results, input_fasta_path)
            
            # Return in the format expected by AlphaFold's run_msa_tool
            return {'a3m': a3m_content}, ''
    
    def _count_hits(self, result_file: str) -> int:
        """Count number of hits in result file."""
        try:
            with open(result_file, 'r') as f:
                return sum(1 for line in f if line.strip())
        except:
            return 0
    
    def _convert_to_a3m(self, result_files: Sequence[str], query_fasta: str) -> str:
        """Convert MMseqs2 results to A3M format."""
        # Read query sequence
        with open(query_fasta, 'r') as f:
            lines = f.readlines()
            query_header = next((line.strip() for line in lines if line.startswith('>')), '>query')
            query_seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
        
        sequences = []
        sequences.append((query_header, query_seq))  # Add query sequence first
        seen_sequences = {query_seq}
        
        # Collect all hits for sorting
        all_hits = []
        
        for result_file in result_files:
            if not os.path.exists(result_file):
                continue
                
            with open(result_file, 'r') as f:
                for line in f:
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 14:
                            try:
                                query_aln = parts[12]   # Aligned query sequence
                                target_aln = parts[13]  # Aligned target sequence
                                target_id = parts[1]
                                evalue = float(parts[10])
                                
                                # Get ungapped sequence for deduplication
                                target_seq_ungapped = target_aln.replace('-', '').upper()
                                
                                if evalue <= 1e-3 and target_seq_ungapped not in seen_sequences:
                                    all_hits.append((evalue, target_id, query_aln, target_aln, target_seq_ungapped))
                            except (ValueError, IndexError) as e:
                                continue
        
        # Sort by E-value (best first)
        all_hits.sort(key=lambda x: x[0])
        
        # Add hits up to max_sequences
        for evalue, target_id, query_aln, target_aln, target_seq_ungapped in all_hits[:self.max_sequences - 1]:
            # Convert to A3M format
            a3m_seq = self._convert_to_a3m_format(query_aln, target_aln)
            if a3m_seq:
                sequences.append((f'>{target_id}', a3m_seq))
                seen_sequences.add(target_seq_ungapped)
        
        # Format as A3M with proper line breaks
        a3m_lines = []
        for header, seq in sequences:
            a3m_lines.append(header)
            # Break long sequences into 80-character lines
            for i in range(0, len(seq), 80):
                a3m_lines.append(seq[i:i+80])
        
        a3m_content = '\n'.join(a3m_lines)
        print(f"Generated A3M with {len(sequences)} unique sequences")
        return a3m_content

    def _convert_to_a3m_format(self, query_aln: str, target_aln: str) -> str:
        """Convert aligned sequences to A3M format."""
        a3m_seq = []
        for q, t in zip(query_aln, target_aln):
            if q == '-':
                # Insertion in target relative to query - lowercase
                if t != '-':
                    a3m_seq.append(t.lower())
            elif t == '-':
                # Deletion in target - skip
                continue
            else:
                # Match/mismatch - uppercase
                a3m_seq.append(t.upper())
        return ''.join(a3m_seq)


# For compatibility, create an HHBlits-like interface
class HHBlits:
    """MMseqs2-based replacement for HHBlits with same interface."""
    
    def __init__(self, binary_path: str, databases: Sequence[str], n_cpu: int = 8):
        self.mmseqs_runner = Mmseqs2HhblitsReplacement(
            binary_path=binary_path,
            databases=databases,
            n_cpu=n_cpu,
            use_gpu=True  # Set to false if you don't have GPU-enabled MMseqs2 :)
        )
    
    def query(self, input_fasta_path: str) -> Mapping[str, Any]:
        return self.mmseqs_runner.query(input_fasta_path)