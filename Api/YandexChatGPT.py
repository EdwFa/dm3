from transformers import AutoTokenizer
from auto_gptq import AutoGPTQForCausalLM
from .configs import check_device, torch


model_name = 'fffrrt/ruGPT-3.5-13B-GPTQ'
model_basename = 'gptq_model-4bit-128g'

device = "cpu"
if check_device():
    device = torch.cuda.current_device()

tokenizer = AutoTokenizer.from_pretrained(model_name, use_fast=True)
model = AutoGPTQForCausalLM.from_quantized(model_name,
        model_basename = model_basename,
        use_safetensors=True,
        trust_remote_code=True,
        device=device,
        use_triton=False,
        quantize_config=None)


def getResponse(request):
    request = f'{request} \n1'
    encoded_input = tokenizer(request, return_tensors='pt').to(device)
    output = model.generate(
        **encoded_input,
        num_beams=4,
        max_new_tokens=160,
        no_repeat_ngram_size=2,
        # num_return_sequences=5,
        # do_sample=True
    )

    print(tokenizer.decode(output[0], skip_special_tokens=True))
    return tokenizer.decode(output[0], skip_special_tokens=True)