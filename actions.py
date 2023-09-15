def get_actions():
	return ACTIONS


def download(context, args):
	context.instrument.download(context, args)


def run_pipeline(context, args):
	context.instrument.run_pipeline(context, args)



ACTIONS = {
	'download': download,
	'pipeline': run_pipeline,
}
