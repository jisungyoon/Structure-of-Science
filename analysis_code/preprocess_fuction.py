def refine_ll_title(ll_decoded, prefix):
    revised_string = ''
    if ll_decoded.startswith(prefix):
        revised_string = ll_decoded[len(prefix):].replace(' ','_')
    else:
        revised_string = ll_decoded.replace(' ','_')
        
    return revised_string