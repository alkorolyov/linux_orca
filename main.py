import asyncio
import time
import logging

logging.basicConfig(
    format='[%(asctime)s] [%(levelname)s] %(message)s',
    datefmt='%d-%m-%Y %I:%M:%S',
    encoding='utf-8',
    level=logging.INFO)


async def run_script(script, script_id, semaphore):
    try:
        async with semaphore:
            process = await asyncio.create_subprocess_shell(
                script, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE
            )
            stdout, stderr = await process.communicate()
            return_code = process.returncode

            # print(f"Script {script_id} completed with return code {return_code}")
            return return_code
    except Exception as e:
        # print(f"Script {script_id} failed with error: {e}")
        return 1

async def run_scripts(script_queue, queue_size):
    semaphore = asyncio.Semaphore(queue_size)

    # Start scripts
    tasks = []
    finished_scripts = 0
    failed_scripts = 0
    while True:
        # Start new scripts if there is room
        while len(tasks) < queue_size and not script_queue.empty():
            script, script_id = await script_queue.get()
            task = asyncio.create_task(run_script(script, script_id, semaphore))
            tasks.append(task)

        # Check for finished tasks
        for task in tasks:
            if task.done():
                if task.result() == 0:
                    finished_scripts += 1
                else:
                    failed_scripts += 1
                tasks.remove(task)

        # Print status update with timestamp
        logging.info(f"Running: {len(tasks)}, Finished: {finished_scripts}, Failed: {failed_scripts}, Waiting: {script_queue.qsize()}")
        await asyncio.sleep(1)

        # # Break the loop if there are no more tasks and the script queue is empty
        if len(tasks) == 0 and script_queue.empty():
            break

    # Wait for all tasks to complete before exiting
    await asyncio.gather(*tasks)

async def main(scripts):
    # Create a script queue with unique script IDs
    script_queue = asyncio.Queue()
    for script_id, script in enumerate(scripts):
        script_queue.put_nowait((script, script_id))

    # Run scripts with a maximum of 3 running concurrently
    await run_scripts(script_queue, queue_size=3)

if __name__ == "__main__":
    # Example scripts
    scripts = [
        "echo 'Script 0'; sleep 1",
        "echo 'Script 1'; sleep 10",
        "echo 'Script 2'; sleep 2",
        "echo 'Script 3'; sleep 4",
        "echo 'Script 3'; sleep 2",
        # Add more scripts as needed
    ]

    asyncio.run(main(scripts))