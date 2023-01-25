
class ChannelInfo:
    def __init__(self, channel, subarray, resolution, min_wave, max_wave):
        self.channel = channel
        self.subarray = subarray
        self.resolution = resolution
        self.min_wave = min_wave
        self.max_wave = max_wave

resolutions = {
    1 : {
            'short':  3630.0,
            'medium': 3620.0,
            'long': 3590.0,
        },
    2 : {
            'short': 3520.0,
            'medium': 3290.0,
            'long':  3250.0,
        },
    3 : {
            'short': 3010.0,
            'medium': 2370.0,
            'long': 2470.0,
        },
    4 : {
            'short': 1680.0,
            'medium': 1695.0,
            'long': 1430.0,
        },
}


class LineInfo:
    def __init__(self, name, wave, channel, subarray):
        self.name = name
        self.wave = wave
        self.channel = channel
        self.subarray = subarray

lines = [
    LineInfo('MgIV_4', 44871., 1, 'short'),
    LineInfo('ArVI_4', 45280., 1, 'short'),
    LineInfo('FeII_4', 48890., 1, 'short'),
    LineInfo('H200S8_5', 50530., 1, 'short'),
    LineInfo('FeII_5', 53400., 1, 'short'),
    LineInfo('FeVIII_5', 54470., 1, 'short'),
    LineInfo('MgVII_5', 55030., 1, 'short'),
    LineInfo('H200S7_5', 55110., 1, 'short'),
    LineInfo('MgV_5', 56100., 1, 'medium'),
    LineInfo('H200S6_6', 61090., 1, 'medium'),
    LineInfo('PAH62_6', 62000., 1, 'medium'),
    LineInfo('ArIII_6', 63670., 1, 'medium'),
    LineInfo('NiII_6', 66360., 1, 'long'),
    LineInfo('FeII_6', 67210., 1, 'long'),
    LineInfo('H200S5_6', 69090., 1, 'long'),
    LineInfo('ArII_6', 69850., 1, 'long'),
    LineInfo('NaIII_7', 73180., 1, 'long'),
    LineInfo('HI63_7', 74600., 2, 'short'),
    LineInfo('NeVI_7', 76520., 2, 'short'),
    LineInfo('PAH77_7', 77000., 2, 'short'),
    LineInfo('FeVII_7', 78140., 2, 'short'),
    LineInfo('ArV_7', 79020., 2, 'short'),
    LineInfo('H200S4_8', 80260., 2, 'short'),
    LineInfo('ArIII_8', 89910., 2, 'medium'),
    LineInfo('MgVII_9', 90090., 2, 'medium'),
    LineInfo('FeVII_9', 95270., 2, 'medium'),
    LineInfo('H200S3_9', 96650., 2, 'medium'),
    LineInfo('SIV_10', 105100., 2, 'long'),
    LineInfo('PAH113_11', 113000., 2, 'long'),
    LineInfo('ClIV_11', 117630., 3, 'short'),
    LineInfo('SIII_12', 120000., 3, 'short'),
    LineInfo('H200S2_12', 122800., 3, 'short'),
    LineInfo('HI76_12', 123700., 3, 'short'),
    LineInfo('NeII_12', 128100., 3, 'short'),
    LineInfo('ArV_13', 131000., 3, 'short'),
    LineInfo('MgV_13', 135200., 3, 'medium'),
    LineInfo('NeV_14', 143200., 3, 'medium'),
    LineInfo('ClII_14', 143700., 3, 'medium'),
    LineInfo('NeIII_15', 155600., 3, 'medium'),
    LineInfo('H200S1_17', 170300., 3, 'long'),
    LineInfo('FeII_17', 179400., 4, 'short'),
    LineInfo('SIII_18', 187100., 4, 'short'),
    LineInfo('FeVI_19', 195527., 4, 'short'),
    LineInfo('ClIV_20', 203200., 4, 'short'), # correct subarray?
    LineInfo('ArIII_21', 218250., 4, 'short'), # correct subarray?
    LineInfo('FeIII_22', 229250., 4, 'short'), # correct subarray?
    LineInfo('NeV_24', 243200., 4, 'long'),
    LineInfo('OIV_25', 258900., 4, 'long'),
    LineInfo('FeII_25', 259900., 4, 'long'),
]

