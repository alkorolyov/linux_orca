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


def log_output(stdout, stderr):
    if stdout.decode():
        logging.debug(f"STDOUT: {stdout.decode()}")
    if stderr.decode():
        logging.debug(f"STDERR: {stderr.decode()}")


async def run_job(job, semaphore):
    cmd, job_id = job['cmd'], job['job_id']
    try:
        async with semaphore:
            process = await asyncio.create_subprocess_shell(
                cmd, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE
            )
            stdout, stderr = await process.communicate()
            return_code = process.returncode

            logging.debug(f"Job {job_id} completed with return code {return_code}")
            log_output(stdout, stderr)
            return return_code

    except Exception as e:
        logging.warning(f"Job {job_id} failed with error: {e}")
        log_output(stdout, stderr)
        return 1


async def main_loop(job_queue, max_parallel):
    semaphore = asyncio.Semaphore(max_parallel)

    # Start jobs
    tasks = []
    finished_jobs = 0
    failed_jobs = 0
    start_time = time.time()
    while True:
        # Start new jobs if there is room
        while len(tasks) < max_parallel and not job_queue.empty():
            job = await job_queue.get()
            task = asyncio.create_task(run_job(job, semaphore))
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
        # logging.info(f"Running: {len(tasks)}, Finished: {finished_jobs}, Failed: {failed_jobs}, Waiting: {job_queue.qsize()}")
        print(" " * 120, end='')
        print(f"\r[ {time.time() - start_time:2.2f}s ] Running: {len(tasks)}, Finished: {finished_jobs}, Failed: {failed_jobs}, Waiting: {job_queue.qsize()}", end='')
        await asyncio.sleep(1)

        # # Break the loop if there is no more tasks and the job queue is empty
        if not tasks and job_queue.empty():
            break

    # Wait for all tasks to complete before exiting
    await asyncio.gather(*tasks)


def main(commands, max_parallel):
    # Create a job queue with unique job IDs
    job_queue = asyncio.Queue()
    for job_id, cmd in enumerate(commands):
        job = {'cmd': cmd, 'job_id': job_id}
        job_queue.put_nowait(job)
    # Run jobs with a maximum `max_parallel` in parallel
    asyncio.run(main_loop(job_queue, max_parallel=max_parallel))

if __name__ == "__main__":

    logging.getLogger().setLevel(logging.INFO)

    # num_jobs = 24
    # cores_per_job = 4
    #
    # commands = [f"python job.py {cores_per_job}"] * num_jobs
    # max_jobs = num_jobs // cores_per_job
    # main(commands, max_parallel=max_jobs)

    num_jobs = 24
    for cores_per_job in [1, 2, 3, 4]:
        max_jobs = num_jobs // cores_per_job
        print('\ncores_per_job:', cores_per_job, 'max_jobs:', max_jobs)
        commands = [f"python job.py {cores_per_job}"] * num_jobs
        main(commands, max_parallel=max_jobs)

    # Example list of bash commands

    # commands = [
    #     # "echo 'job 0'; sleep 1",
    #     # "echo 'job 1'; sleep 5",
    #     # "echo 'job 2'; sleep 2; exit 1",
    #     # "echo 'job 3'; sleep 3",
    #     # "echo 'job 4'; sleep 1",
    #     # Add more jobs as needed
    # ]
