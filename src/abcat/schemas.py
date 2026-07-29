from enum import Enum

from pydantic import BaseModel, Field


class ChainType(str, Enum):
    HEAVY = "Heavy"
    KAPPA = "Kappa"
    LAMBDA = "Lambda"
    UNKNOWN = "Unknown"


class CDRRegion(BaseModel):
    name: str  # "CDR1", "CDR2", "CDR3"
    sequence: str
    start_pos: str  # IMGT start coordinate (e.g. "27")
    end_pos: str  # IMGT end coordinate (e.g. "38")
    length: int


class FrameworkRegion(BaseModel):
    name: str  # "FR1", "FR2", "FR3", "FR4"
    sequence: str


class VAnalysisResult(BaseModel):
    chain_type: ChainType
    species: str | None = "human"
    scheme: str = "imgt"
    v_gene_allele: str | None = None
    j_gene_allele: str | None = None
    v_domain_sequence: str
    numbering: dict[str, str] = Field(
        default_factory=dict
    )  # {"27": "G", "28": "F", ...}
    cdrs: dict[str, CDRRegion] = Field(default_factory=dict)
    frameworks: dict[str, FrameworkRegion] = Field(default_factory=dict)
    confidence: float = 1.0
    error_message: str | None = None


class AllotypeMarkerCall(BaseModel):
    allotype: str
    opposing_allotype: str | None = None
    category: str = "allotype"  # "allotype" or "isoallotype"
    domain: str  # "CH1", "CH2", "CH3", "CL"
    status: str  # "present", "absent", "inconclusive"
    matched_residues: dict[str, str] = Field(default_factory=dict)  # {"EU_214": "K"}


class CAnalysisResult(BaseModel):
    c_region_sequence: str
    isotype: str  # "IgG", "IgA", "IgM", "IgE", "IgD", "Kappa", "Lambda"
    subclass: (
        str  # "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2", "IGKC", "IGLC1", ...
    )
    matched_c_gene: str  # "IGHG1", "IGKC", etc.
    alignment_identity: float
    allotypes: list[AllotypeMarkerCall] = Field(default_factory=list)
    isoallotypes: list[AllotypeMarkerCall] = Field(default_factory=list)
    allotype_summary: str = "None"


class FullChainAnalysis(BaseModel):
    sequence_id: str
    full_sequence: str
    v_analysis: VAnalysisResult | None = None
    c_analysis: CAnalysisResult | None = None
    success: bool = True
    warnings: list[str] = Field(default_factory=list)


class BatchAnalysisResult(BaseModel):
    total_sequences: int
    successful_analyses: int
    results: list[FullChainAnalysis]
