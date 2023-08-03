import torch


def check_device():
    if not torch.cuda.is_available():
        return False

    print(torch.cuda.device_count())
    print(torch.cuda.current_device())
    return True