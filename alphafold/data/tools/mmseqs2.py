"""MMseqs2 wrapper to replace HHblits in AlphaFold MSA generation."""

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
                 max_sequences: int = 10000):
        """Initialize MMseqs2 runner."""
        self.binary_path = binary_path
        self.databases = databases
        self.n_cpu = n_cpu
        self.sensitivity = sensitivity
        self.max_sequences = max_sequences
        
        # Verify binary exists
        if not os.path.exists(self.binary_path):
            raise FileNotFoundError(f"MMseqs2 binary not found at {self.binary_path}")
        
        # Verify databases exist
        for db in self.databases:
            if not os.path.exists(db):
                print(f"Warning: Database not found at {db}")
    
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
                    '--format-output', 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qseq,tseq',
                    '-e', '1e-3',  # E-value threshold
                    '--comp-bias-corr', '1',
                    '--remove-tmp-files'  # Clean up temporary files
                ]
                
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
            return {'a3m': a3m_content}
    
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
            query_seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
            query_header = next((line.strip() for line in lines if line.startswith('>')), '>query')
        
        sequences = []
        sequences.append((query_header, query_seq))  # Add query sequence first
        seen_sequences = {query_seq}  # Track unique sequences
        
        # Parse MMseqs2 results
        for result_file in result_files:
            if not os.path.exists(result_file):
                continue
                
            with open(result_file, 'r') as f:
                for line in f:
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 14:
                            target_seq = parts[13]  # Target sequence with gaps
                            target_id = parts[1]    # Target identifier
                            evalue = float(parts[10])  # E-value
                            
                            if target_seq and target_seq not in seen_sequences and evalue <= 1e-3:
                                # Convert dots to gaps for A3M format
                                a3m_seq = target_seq.replace('.', '-')
                                sequences.append((f'>{target_id}', a3m_seq))
                                seen_sequences.add(target_seq)
        
        # Generate A3M format string
        a3m_lines = []
        for header, seq in sequences:
            a3m_lines.append(header)
            a3m_lines.append(seq)
        
        a3m_content = '\n'.join(a3m_lines)
        print(f"Generated A3M with {len(sequences)} unique sequences")
        return a3m_content


# For compatibility, create an HHBlits-like interface
class HHBlits:
    """MMseqs2-based replacement for HHBlits with same interface."""
    
    def __init__(self, binary_path: str, databases: Sequence[str], n_cpu: int = 8):
        self.mmseqs_runner = Mmseqs2HhblitsReplacement(
            binary_path=binary_path,
            databases=databases,
            n_cpu=n_cpu
        )
    
    def query(self, input_fasta_path: str) -> Mapping[str, Any]:
        return self.mmseqs_runner.query(input_fasta_path)