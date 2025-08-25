import os
import sys
import time
import logging
from rich.console import Console
from datetime import datetime
class Logger:
   def __init__(self, log_file):
       self.console = Console()
       self.log_file = log_file
       self.logger = logging.getLogger('Logger')
       self.logger.setLevel(logging.INFO)
       handler = logging.FileHandler(self.log_file)
       handler.setLevel(logging.INFO)
       formatter = logging.Formatter('%(asctime)s - %(message)s')
       handler.setFormatter(formatter)
       self.logger.addHandler(handler)

   def log(self, message):
       self.console.log(message)
       self.logger.info(message)

   def start_timer(self):
       self.start_time = time.time()

   def stop_timer(self):
       end_time = time.time()
       execution_time = end_time - self.start_time
       self.log(f"Execution time: {execution_time} seconds")


logger = Logger('vir_disc_pipeline.log')