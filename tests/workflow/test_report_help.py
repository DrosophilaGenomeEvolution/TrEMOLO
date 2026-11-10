import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest
from tests.workflow import test_build_quarto_report as fixtures

ROOT = Path(__file__).resolve().parents[2]


class ReportHelpTests(unittest.TestCase):
    def test_cli_writes_offline_help_and_context_links_resolve(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            infos,index,manifest,_,_=fixtures.BuildQuartoReportTests().fixture(root)
            out=root/'REPORT'
            command=[sys.executable,str(ROOT/'lib/python/workflow/build_quarto_report.py')]
            for key,value in {'te-infos':infos,'genome-index':index,'input-manifest':manifest,
                              'template':ROOT/'report/quarto/report.qmd','style':ROOT/'report/quarto/report.css',
                              'script':ROOT/'report/quarto/report.js','output-qmd':out/'report.qmd','output-json':out/'report-data.json'}.items():
                command.extend(['--'+key,str(value)])
            result=subprocess.run(command,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            help_page=(out/'report-help.html').read_text()
            anchors=set(re.findall(r'id="([^"]+)"',help_page))
            for anchor in re.findall(r'report-help\.html#([a-z-]+)',(out/'report.qmd').read_text()):
                self.assertIn(anchor,anchors)
            for term in ['nested_candidate','same_family_overlap','accepted','partial','full_length']:
                self.assertIn('<code>'+term+'</code>',help_page)
            self.assertNotRegex(help_page,r'(?:src|href)="https?://')
            data=json.loads((out/'report-data.json').read_text())
            self.assertFalse(data['sv_counts_available'])
