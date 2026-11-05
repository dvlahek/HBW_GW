#!/usr/bin/env python3
"""Frozen local C8 Taylor surrogates from genuine reconstructed q=.3 C1-C7 nodes.

Compressed storage is byte-for-byte reversible to JSON. SHA256 of compressed
payload is fixed. All underlying q=.3 NPZ SHA256 and polynomial errors are
in the decoded JSON, with separate local replication sources. The original
new-q NPZ is NOT committed, so this is exploratory analytical testing only.
"""
import base64,hashlib,json,zlib
DATA_B64="""eNrtmU1vG0cShv8Lz/Sgu6r6y9dcggAL5LiAYRAUOVoTUCwvJTv2Gv7v+1TPkCYlJhCze9gsEkOGQ3bPdFe99dbTra+L99niah0Wr78u/rl4HQZdLj68+/Kw26zvVg8fxs3j3v/xbi0pL14vmoxZW9FtiVbsxlqstxvZSEtNpKTcVOw2t+1t2ba63m7X21i13mgbiwbRtlgu3n/eLV7z0uVi/X7z7n6/ur+9fRgfV/tx9bf7X8Z/rH0VYbl4GMftanP/y4e78fPxmzcssLQSouaaRWMzXb5iOI+uySrf1axV3i4Xm/XHzbsvq81uv7kbV/v1dvfxoT851OXicf3lzt+83477xesovO1xv9s8rj6t73bb3eOXswnG0+7H29vdZje+f3zwUP38008s5o1/G1JokZc2bbmFVNvSP81VrEmKWgt/cnm7fPOqj9aSVRgeY2stJcvLV3EIicD5Z0FFYszjq5CYYQPbEZ4aW4otKx/npQ7RNDcxSWal5TgPlkFD0tYkWq2W1KwP5+lK7HMIxEpTSY2P+3JsKJqziPETm2gufFOXr2QgZWYiqlGKTp/3GWXIUZrUEkpha2WeoANht1Yk+fNT8hU1Hx+HGEvyDWgONRWVPiEOrTRLibCgG8u+3cp4HUxKjspfIZv4QtsyscoUCAyh1VL7I/zhPERrjQSAz2vVoH04qzfeI7WiCVYT6vgqxh6eYiyjemyqxkJ0YlgShIqeGsmQUoPF4h/72usgBFGSqPX9ij+nR7Oy9ZgLY3OJ6TDeBtdAkZCTRvVpfUIZLFaSXiIRMvXHqw/PAwFg/4iYkinGZpEhASNY2cxjVVF5f+tb5PzzT38/CK76902thZLJh3a9WclMak15IAn6rjcR4uzVyeatpjwNDxRx8WS5CFl06eHnI95eooRYqnVh+Y6l1ZZyzoHlmtlBnYQI8UtOxcgLofXUlEk/hCagUaHouxyyj1dCJJ5iyRnjwC76hDY0yS5pIoXMVWc5MJ7/tUDypRRmTRluQ/EkVhRd0fpxuAxJDEcIBJZtREuz2Iz4Zi1oPLLWvhwfX9hXYu2xsVx0Uif92IBsYhGLhMyyZ2wSsyIIAiuRT0NWr0YXBHrL1BrrzK7/XrttyrAQfcqanFCoWSbBIfLi1Y5bEs0a2kFAcSipJsbmGhI/uUuL9QerFqh1MkaxVJ0ntKGiWiq7Ejhe0qVFuWBEiacjXnw6zeonChijoIPQ685rmtXXIdaW8adAAfGa7CrscvvxoDaLJjzeSGXGjVK3W0KJ2HgL9eERP6gtsaQYSvAluE/3wQG7wxqpO3QVifdhdEQMJNEQibnrdWVSWdQOMSabHlgGd8ESQkWAElPwmXWScYoBA2AogkXT9TAaGWVFA7yNZ5PhaSVIicLCrUJlBzV6p3hDCjEn9/DYSEEOWifl56EFfBFrSVhATekgfKyWzSTkV/i0pu49iVTV6gpz6ycIxOagfJIYvTxxMk9XjZOV46cBp3LXzgkBzcN5PGWFgVEqCTn3AiJX5jWlBMa/bnV2ZQTr7iU4gsauoKkMWT1m3chsa9Zaf2WZFkObwr1rJEzegua6ZfkJVcbi+eO7vt/Z9yVmPuS/FFlnL7kl/QPlYRU0EaKT21xYbFZ91d4raMk1Tq4vQ5aGh2CFfIFzTYX79ttyMe739/vvffXr4uH+434zrn4oK/Djw/iv+9129elhNX5e06N/qKsNfXjcgwx3dO8B68CejL+7N07hmviBfv5x/D71w27cjL/uHsYXTAUNdp/Wj7tPvz0/Ye4RZ8cDsb8yNdzFr+zl8Ugu78a77f3HR3/G97302XQYQXU0WjRfunNE+zZb/bUxoApFWAcWKvzkK2JwnOpIFWZpXhEDOoMWXBGsCfw7XRcFtIWgCIJUTKu2dojCj9cGQclmTNENwuYqjuFFMTjMDPhppB+3dDL1BTHIA1WNkzpGBYnT5JcGwAEHH0zeJFKU1C372zdC8B4CeDmal7Su67YOa1uH8WZLY8QQbtalrAswqGvsUNbbti43N+lmHbdxY7eMqkT+Rra8ZkZzXvp/g+ZWEDcmTVgxYBrejOZSXaoC3sbUTtCczghqQJXkU8WbO3YX6ROVCEJY5QjmMAX8QkNlXxNqdDKHT3kIUWBW7QQ+oTloAJzAWdH5p03jKT11UHUPtmyhyGy4gA/mGTyKMCVoYge21eCI1PyYYNjOwaDTkGDhmAHKguf33gCYV55DQ8LPOYvYjMg+vuHFwcQ9i0hIf3FzVPKI4d4K8fph5IhWVAfLYV9OQGFajg3e+4gC7sWSih7J3NmeUw9ogb5CPyb0pzcHXXQYaD21A3KYwc1ZlpfXQmYsHUhJnaVJH8DIDD1h7ZS8SjEtNKW9XDubEwI/yGjiRSG34wSgi1CjxMYbzI9I04w0OOnxXuCY53SalxneCIFQ1fRUTl4djLxfG+DjMOBYV8sRl87o3HslfRJyQHUTnRMF2gRdF55q0wlgonPIiFbsgoMPEN2ENUBvoTOAZcgpmaYuiupbYwOcnQhMO+I5SQNZc3DKBwMOlMLRBo1IdpwCbpod8JyXAW7msQLJ8gE7+Ib+DFGzzKROOH0CsQMkQUbzGTmmA563gRcm9opCnWomPHdWcKBWCo3xSY9UALxxiENETkLSjxGO5+IVZ75n4NGyzOOxxZY7PZB6DpHpUAJkxS0FyOrVakfV1Y51cDJv9gPUpCI+d9VF13lyRjwCOkXMMcSLVaPleORzvLj4SzjIAvSdiWdAJwAEH9XR8KDe6Q0slBV6DXDIdlHOyI2GmI71kkiD7zwxaAhc89M4Ry8eQKrLPJxjDesA2dSTzoFoVhwh4ZRekQ9ngN4dLwK6BwPzIeGnhE7eOQtx2s3pHNEdJDEeQOwU0fEhPz4hez1ndMATWRROH/U7pKvXJTMxznNIB7gJW/MjzAmlEwO/a/ATzCmli/dc3BRb9LuQE0rn8aSOU8xk0TOlG2lUAFY64xwpHT0rKiUCTyCdXJF680P2KaQrzge+1pLtCaSjLvLGaXT6ZoZ0TNQD5rV4RunFfYtQkwY9wXQONS5dUnmO6d5LuqUXO6V0MufnZ8KUnlA6he3uRAPTU0r3klQPTsxyTul+mkAJbeoOM6X38CMsVpbsHNMxD3oYGo31FNM9OtRMjiH/dzC9Up/gLoXP8ScfDh0vobPTmU5KpzNfAGfu/ULMWiQQ"""
COMPRESSED_SHA256="08bca5352c795ddf08deb410d7e4f62380d23d28f3e86cbb9cdfde0b5d8b3ea0"
def cases():
    raw=base64.b64decode(DATA_B64)
    if hashlib.sha256(raw).hexdigest()!=COMPRESSED_SHA256:
        raise ValueError("Spectral Taylor provenance checksum mismatch: actual="+hashlib.sha256(raw).hexdigest()+", expected="+COMPRESSED_SHA256+", bytes="+str(len(raw)))
    data=json.loads(zlib.decompress(raw).decode())
    if set(data)!={"n641_a0","n801_a0","n801_am0005"}:
        raise ValueError("Frozen three-case physical-source grid altered")
    for name,x in data.items():
        if x["q"]!=.3 or x["taylor_order"]!=12 or x["strict_validity_radius"]!=.004:
            raise ValueError("Unexpected surrogate model/grid")
        if x["nxi"] not in (641,801) or len(x["physical_spectral_sha256"])!=64:
            raise ValueError("Original q0.3 physical source provenance missing")
        for v in x["errors"].values():
            if (v["seed_value_vs_exact_piecewise_rel"]>1e-4 or
                v["seed_derivative_vs_exact_piecewise_rel"]>1e-4 or
                v["worst_complex_heldout_vs_trapezoid_rel"]>1e-7):
                raise ValueError("Taylor surrogate accuracy has not passed")
    return data
