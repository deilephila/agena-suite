import agena_panels_MFE #скачиваем весь модуль

from fastapi import FastAPI, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi import Request
import json
from pydantic import BaseModel
from enum import Enum
from typing import Optional

app = FastAPI()

origins = [
    "http://127.0.0.1",
    "http://127.0.0.1:8000"
    "http://127.0.0.1:8080",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class classData(BaseModel):
    UEP: dict
    Panel: dict
    TM: dict

@app.post("/addUEP")
async def create_data(predata: classData): 
    dataUEP = predata.UEP
    dataPanel = predata.Panel
    dataTM = predata.TM
    afterMFE = agena_panels_MFE.addUEP(dataUEP,dataPanel,dataTM)
    return afterMFE


app.mount("/", StaticFiles(directory="static",html = True), name="static")