import pathlib
from pydantic import AfterValidator, BaseModel, Field, PlainSerializer
from typing import Annotated, ClassVar, Literal, Type, types


from daedalus.instruments.instrument_common import Instrument, SUPPORTED_INSTRUMENTS
from daedalus.instruments.nirspec_ifu import NIRSpec_IFU
from daedalus.instruments.miri_ifu import MIRI_IFU
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

    instrument_class: ClassVar[Type[Instrument]] = Instrument

    def create_instrument(self, config):
        return self.instrument_class(config)


class NIRSpecIFUInstrumentConfig(InstrumentConfig):
    name: Literal['nirspec_ifu'] = 'nirspec_ifu'
    obs: int = Field(description='The observation number for the target found in the program PDF.')

    instrument_class: ClassVar[Type[Instrument]] = NIRSpec_IFU


class MIRIInstrumentConfig(InstrumentConfig):
    name: Literal['miri'] = 'miri'
    sci_obs: int = Field(description='The science observation number for the target found in the program PDF.')
    bkgd_obs: int = Field(description='The background observation number for the target found in the program PDF.')

    instrument_class: ClassVar[Type[Instrument]] = MIRI_IFU


class TargetConfig(BaseModel):
    name: str = Field(description='Name of target. This does not have to be any particularly referenced name, it is for your benefit.')
    program: int = Field(description='The program ID of the target')
    instrument: NIRSpecIFUInstrumentConfig | MIRIInstrumentConfig = Field(discriminator='name', description='Configuration of the instrument and observation of the target for this build.')


class ActionConfig(BaseModel):
    name: Literal[*list(Instrument.ACTIONS.keys())] = Field(description='The name of the action')
    overwrite: bool = Field(False, description='Overwrite existing files.')


class DownloadConfig(ActionConfig):
    name: Literal['download'] = 'download'
    mast_token: str | None = Field(None, description='MAST API token. Only needed to download proprietary data that you have access to.')


class PipelineStageOptions(BaseModel):
    stage: Literal[1,2,3] = 1
    steps: dict = Field(default_factory=dict, description='Additional options to send to specific steps in this stage.')
    skip: bool = Field(False, description='Skip this whole stage')


class PipelineStage1Options(PipelineStageOptions):
    flag_snowballs: bool = Field(True, description='TODO') # JWST 1.20.2 / CRDS 13.1.10 default
    clean_noise: bool = Field(False, description='TODO') # JWST 1.20.2 / CRDS 13.1.10 default


class PipelineStage2Options(PipelineStageOptions):
    stage: Literal[2] = 2
    cube_build: bool = Field(False, description='TODO') # 


class PipelineStage3Options(PipelineStageOptions):
    stage: Literal[3] = 3
    outlier_detection: bool = Field(False, description='TODO')


class PipelineConfig(ActionConfig):
    name: Literal['pipeline'] = 'pipeline'
    multiprocess: bool = Field(False, description='Switch to run the pipeline build in multiprocess mode.')
    nprocesses: int = Field(1, description='The number of processes to use when running in multiprocess mode.')
    stage1: PipelineStage1Options = Field(PipelineStage1Options(), description='Options for Stage 1 of the pipeline build.')
    stage2: PipelineStage2Options = Field(PipelineStage2Options(), description='Options for Stage 2 of the pipeline build.')
    stage3: PipelineStage3Options = Field(PipelineStage3Options(), description='Options for Stage 3 of the pipeline build.')


class Config(BaseModel):
    product_name: str = Field(description='Unique product name for this build. It can be anything, but will be used to refer to this build.')
    target: TargetConfig = Field(description='Configuration for observed target of this build.')
    actions: list[Annotated[DownloadConfig | PipelineConfig, Field(discriminator='name')]] = Field(default_factory=list, description='The actions to be performed during this build.')
    uncal_dir: Path = Field(None, description='Directory of uncalibrated data. If one is not provided, Daedalus will create one for you and download the necessary files.')
    stage2_dir: Path = Field(None, description='Directory of data processed by stage 1, ready for stage 2. If one is not provided, Daedalus will create one for you and download or build the necessary files.')
    stage3_dir: Path = Field(None, description='Directory of data processed by stage 2, ready for stage 3. If one is not provided, Daedalus will create one for you and download or build the necessary files.')
    output_dir: Path = Field(description='Daedalus build output directory. Note that the final output directory will be: `<output_dir>/<program_id>/<target_name>/<instrument_name>/<product_name>/`. Therefore, the `output_dir` passed in the configuration file can be a general Daedalus build directory for all of your targets and build products.')
    crds_cache: Path = Field(None, description='CRDS cache directory. Will default to some directory as determined by the CRDS module. Warning! These directories can get up to 10s of GBs.')
    # TODO: save_intermediate


def generate_docs(model, prefix=''):
    doc_str = ''

    for field_name, field in model.model_fields.items():

        doc_str += '## `%s%s`\n'%(prefix,field_name)

        def get_type(arg):
            if hasattr(arg, '__name__'): return arg.__name__
            return str(arg)

        if (isinstance(field.annotation, types.UnionType)) or (field.annotation.__name__ == 'Literal') or (field.annotation.__name__ == 'Union'):
            type_str = ' | '.join([get_type(a) for a in field.annotation.__args__])
        else:
            type_str = get_type(field.annotation)
        doc_str += '*Type:* `%s`<br/>\n'%type_str

        if str(field.default) == 'PydanticUndefined':
            default_str = 'None [Required Field]'
        else:
            default_str = str(field.default)
        doc_str += '*Default:* `%s`<br/>\n'%default_str

        if not field.description is None:
            doc_str += '*Description:* %s\n'%field.description

        if getattr(field.annotation, '__name__', None) == 'Literal':
            doc_str += '*Options:* %s\n'%', '.join(['`%s`'%a for a in field.annotation.__args__])

        doc_str += '\n'

        if getattr(field.annotation, '__name__', None) == 'Literal':
            # Literals get angry about the other conditionals, so just skip 'em
            pass
        elif getattr(field.annotation, '__name__', None) == 'list':
            item_types = field.annotation.__args__[0].__args__[0].__args__ # there's certainly a better way to do this
            for subfield in item_types:
                if issubclass(subfield, BaseModel):
                    doc_str += generate_docs(subfield, prefix=prefix+'%s.'%subfield.__name__)
        elif (isinstance(field.annotation, types.UnionType)) or (field.annotation.__name__ == 'Union'):
            for subfield in field.annotation.__args__:
                if issubclass(subfield, BaseModel):
                    doc_str += generate_docs(subfield, prefix=prefix+'%s.'%subfield.__name__)
        elif issubclass(field.annotation, BaseModel):
            doc_str += generate_docs(field.annotation, prefix=prefix+'%s.'%field_name)


    return doc_str


if __name__ == '__main__':
    doc_dir = pathlib.Path(__file__).resolve().parent.parent.parent.joinpath('docs', 'documentation', 'usage')
    doc_str = generate_docs(Config)
    with open(doc_dir.joinpath('config_options.md'), 'w') as f:
        f.write(doc_str)




