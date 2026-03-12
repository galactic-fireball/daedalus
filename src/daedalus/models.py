import pathlib
from pydantic import AfterValidator, BaseModel, Field, PlainSerializer
from typing import Annotated, Literal


from daedalus.actions import ACTIONS
from daedalus.instruments.instrument_common import SUPPORTED_INSTRUMENTS
from daedalus.utils import REPO_DIR


def to_path(v: pathlib.Path | str) -> pathlib.Path:
	if v is None: return v

	path = pathlib.Path(v)
	if not path.is_absolute():
		path = REPO_DIR.joinpath(path)

	return path


def serialize_path(v: pathlib.Path | str) -> str:
	return str(v)


Path = Annotated[pathlib.Path | str, AfterValidator(to_path), PlainSerializer(serialize_path, return_type=str)]


class InstrumentConfig(BaseModel):
	name: Literal[*SUPPORTED_INSTRUMENTS] = Field(description='The name of the instrument.')


class NIRSpecIFUInstrumentConfig(InstrumentConfig):
	name: Literal['nirspec_ifu'] = 'nirspec_ifu'
	obs: int = Field(description='The observation number for the target found in the program PDF.')


class MIRIInstrumentConfig(InstrumentConfig):
	name: Literal['miri'] = 'miri'
	sci_obs: int = Field(description='The science observation number for the target found in the program PDF.')
	bkgd_obs: int = Field(description='The background observation number for the target found in the program PDF.')


class TargetConfig(BaseModel):
	name: str = Field(description='Name of target.')
	program: int = Field(description='The program ID of the target')
	redshift: float = Field(0.0, description='Redshift of target')
	instrument: NIRSpecIFUInstrumentConfig | MIRIInstrumentConfig = Field(discriminator='name', description='Configuration of the instrument and observation of the target for this build.')


class ActionConfig(BaseModel):
	name: Literal[*ACTIONS] = Field(description='The name of the action')


class DownloadConfig(ActionConfig):
	name: Literal['download'] = 'download'


class PipelineStageOptions(BaseModel):
	stage: Literal[1,2,3] = 1
	steps: dict = Field(default_factory=dict, description='Additional options to send to specific steps in this stage.')


class PipelineStage1Options(PipelineStageOptions):
	flag_snowballs: bool = Field(False, description='TODO')
	clean_rates: bool = Field(False, description='TODO')


class PipelineStage2Options(PipelineStageOptions):
	stage: Literal[2] = 2


class PipelineStage3Options(PipelineStageOptions):
	stage: Literal[3] = 3
	outlier_detection: bool = Field(False, description='TODO')


class PipelineConfig(ActionConfig):
	name: Literal['pipeline'] = 'pipeline'
	multiprocess: bool = Field(False, description='Switch to run the pipeline build in multiprocess mode.')
	stage1: PipelineStage1Options = Field(PipelineStage1Options(), description='Options for Stage 1 of the pipeline build.')
	stage2: PipelineStage2Options = Field(PipelineStage2Options(), description='Options for Stage 2 of the pipeline build.')
	stage3: PipelineStage3Options = Field(PipelineStage3Options(), description='Options for Stage 3 of the pipeline build.')


class Config(BaseModel):
	product_name: str = Field(description='Unique product name for this build.')
	output_dir: Path = Field(description='Directory to put build output.')
	crds_cache: Path = Field(None, description='CRDS cache directory.')
	use_crds_ops: bool = Field(True, description='TODO')
	target: TargetConfig = Field(description='Configuration for observed target of this build.')
	actions: list[Annotated[DownloadConfig | PipelineConfig, Field(discriminator='name')]] = Field(default_factory=list, description='The actions to be performed during this build.')

