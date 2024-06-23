

import os
from ops.argparser import argparser_train

if __name__ == "__main__":
    params = argparser_train()
    gpu_id = params['gpu']
    if gpu_id is not None:
        os.environ["CUDA_VISIBLE_DEVICES"] = gpu_id
    from training.main_worker import main_worker
    main_worker(params)