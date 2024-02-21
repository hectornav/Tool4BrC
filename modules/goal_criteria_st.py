import numpy as np

def goal_nmb(x, pollutant='oa'):
    if pollutant == 'oa':
        text = None if np.abs(x) < 15 else ''
        color = 'gray' if np.abs(x) < 15 else None
        if color == None:
            pass
        else:
            return f"font-style: {text}; background-color: {color};"
            
    elif pollutant == 'bc':
        text = 'bold' if np.abs(x) < 20 else ''
        color = 'gray' if np.abs(x) < 20 else None
        if color == None:
            pass
        else:
            return f"font-style: {text}; background-color: {color};"

def goal_nme(x, pollutant='oa'):
    if pollutant == 'oa':
        text = None if np.abs(x) < 45 else ''
        color = 'gray' if np.abs(x) < 45 else None
        if color == None:
            pass
        else:
            return f"font-style: {text}; background-color: {color};"
    elif pollutant == 'bc':
        text = 'bold' if np.abs(x) < 50 else ''
        color = 'gray' if np.abs(x) < 50 else None
        if color == None:
            pass
        else:
            return f"font-style: {text}; background-color: {color};"

def goal_fb(x, pollutant='oa'):
    if pollutant == 'oa':
        text = None if np.abs(x) < 30 else ''
        color = 'gray' if np.abs(x) < 30 else None
        if color == None:
            pass
        else:
            return f"font-style: {text}; background-color: {color};"
    elif pollutant == 'bc':
        text = 'bold' if np.abs(x) < 30 else ''
        color = 'gray' if np.abs(x) < 30 else None
        if color == None:
            pass
        else:
            return f"font-style: {text}; background-color: {color};"

def criteria_fac2(x, pollutant='oa'):
    if pollutant == 'oa':
        text = None if np.abs(x) >= 50 else ''
        weight= 'bold' if np.abs(x) >= 50 else ''
        return f"font-style: {text}; font-weight:{weight};"
    elif pollutant == 'bc':
        text = 'bold' if np.abs(x) >= 50 else ''
        weight= 'bold' if np.abs(x) >= 50 else ''
        return f"font-style: {text}; font-weight:{weight};"
        
def criteria_nme(x, pollutant='oa'):
    if pollutant == 'oa':
        text = 'bold' if np.abs(x) < 65 else ''
        return f"font-weight: {text}"
    elif pollutant == 'bc':
        text = 'bold' if np.abs(x) < 75 else ''
        return f"font-weight: {text}"

def criteria_nmb(x, pollutant='oa'):
    if pollutant == 'oa':
        text = 'bold' if np.abs(x) < 50 else ''
        return f"font-weight: {text}"
    elif pollutant == 'bc':
        text = 'bold' if np.abs(x) < 40 else ''
        return f"font-weight: {text}"

def criteria_fb(x, pollutant='oa'):
    if pollutant == 'oa':
        text = 'bold' if np.abs(x) < 60 else ''
        return f"font-weight: {text}"
    elif pollutant == 'bc':
        text = 'bold' if np.abs(x) < 60 else ''
        return f"font-weight: {text}"