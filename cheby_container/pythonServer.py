
import time
#from sqs_class import aws_sqs
import logging

#aws_sqs = aws_sqs()

logging.basicConfig(
	format = '%s(asctime)s %(levelname)-8s %(message)s', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S', filename='pythonLog')


while (1):
	'''
	message = aws_sqs.receive_message( )
	if message == None:

		logging.info('No new message in Queue')

	else:
		logging.info("New Message in Queue")
		logging.info(message)

		#for key in message:
		#	print (key)
		if message['MessageAttributes']['Title']['StringValue'] == 'FindOrb':
			print ('Fire UP FindOrb')
		
		#stop

	print ('hello PythonServer is working')
	'''
	time.sleep(60)