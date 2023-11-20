import asyncio
import sys
import time
import logging
import argparse

logging.basicConfig(
    format='[%(asctime)s] [%(levelname)s] %(message)s',
    datefmt='%d-%m-%Y %I:%M:%S',
    encoding='utf-8',
    level=logging.INFO)

async def run_job(cmd, job_id, semaphore):
    try:
        async with semaphore:
            process = await asyncio.create_subprocess_shell(
                cmd, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE
            )
            stdout, stderr = await process.communicate()
            return_code = process.returncode

            logging.debug(f"Job {job_id} completed with return code {return_code}")
            logging.debug(f"STDOUT: {stdout.decode()}")
            return return_code

    except Exception as e:
        logging.warning(f"Job {job_id} failed with error: {e}")
        logging.debug(f"STDOUT: {stdout.decode()}")
        logging.debug(f"STDERR: {stderr.decode()}")
        return 1

async def event_loop(job_queue, queue_size):
    semaphore = asyncio.Semaphore(queue_size)

    # Start jobs
    tasks = []
    finished_jobs = 0
    failed_jobs = 0
    while True:
        # Start new jobs if there is room
        while len(tasks) < queue_size and not job_queue.empty():
            cmd, job_id = await job_queue.get()
            task = asyncio.create_task(run_job(cmd, job_id, semaphore))
            tasks.append(task)

        # Check for finished tasks
        for task in tasks:
            if task.done():
                if task.result() == 0:
                    finished_jobs += 1
                else:
                    failed_jobs += 1
                tasks.remove(task)

        # Print status update with timestamp
        logging.info(f"Running: {len(tasks)}, Finished: {finished_jobs}, Failed: {failed_jobs}, Waiting: {job_queue.qsize()}")
        await asyncio.sleep(1)

        # # Break the loop if there is no more tasks and the job queue is empty
        if not tasks and job_queue.empty():
            break

    # Wait for all tasks to complete before exiting
    await asyncio.gather(*tasks)

def main(commands):
    # Create a job queue with unique job IDs
    job_queue = asyncio.Queue()
    for job_id, cmd in enumerate(commands):
        job_queue.put_nowait((cmd, job_id))

    # Run jobs with a maximum `queue_size` in parallel
    asyncio.run(event_loop(job_queue, queue_size=4))


if __name__ == "__main__":

    # parser = argparse.ArgumentParser(description='Description of your script')
    # # Define command-line arguments
    # parser.add_argument('input_file', help='Path to the input file')
    #
    # # Parse the command-line arguments
    # args = parser.parse_args()
    #
    # # read filepath containing list of commands to run
    # filepath = sys.argv[1]
    # if filepath:
    #     with open(filepath, 'r') as f:
    #         commands = f.readlines()

    # Example list of bash commands
    commands = [
        "echo 'job 0'; sleep 1",
        "echo 'job 1'; sleep 5",
        "echo 'job 2'; sleep 2; exit 1",
        "echo 'job 3'; sleep 3",
        "echo 'job 4'; sleep 1",
        # Add more jobs as needed
    ]

    # asyncio.run(main(commands))
    main(commands)