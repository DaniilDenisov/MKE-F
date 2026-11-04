import unittest
from unittest.mock import patch

from section_solver import Engine, Options, Profile
from section_service.display import export
from section_service.jobs import Jobs


class DisplayLimits(unittest.TestCase):
    def test_individual_limit(self):
        load = [1,0,0,0,0,0]
        k = Engine().solve(Profile("rectangle",2,1),.3,load,Options(levels=1))["kernel"]
        meta,data = export(k,load,limit=1)
        self.assertIsNone(data)
        self.assertEqual(meta["status"],"resource_limit")

    def test_retention_and_expiry(self):
        jobs = Jobs(4)
        with patch("section_service.jobs.RETAINED_LIMIT",10):
            for key in ("old","new"):
                jobs.records[key] = {}
                jobs.finish(key,"completed",{"display":{"status":"ok"}},b"12345678")
            self.assertEqual(jobs.get("old")["result"]["display"]["status"],"evicted")
            self.assertEqual(jobs.field("new"),b"12345678")
            jobs.records["new"]["finished"] -= 901
            jobs.prune()
            self.assertFalse(jobs.fields)
