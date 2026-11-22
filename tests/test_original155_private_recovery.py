"""Input-integrity failures must stop original155 restoration before model use."""
from __future__ import annotations
import hashlib
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest
import zipfile
import zlib
import contextlib
from unittest.mock import patch

MODULE_PATH = Path(__file__).resolve().parents[1] / 'scripts/recover_GEB_v05_original155_private_inputs_20261122.py'


class RecoveryTests(unittest.TestCase):
    def setUp(self):
        self.assertTrue(MODULE_PATH.is_file(), 'Original155 verified recovery implementation is missing')
        spec = importlib.util.spec_from_file_location('recovery155', MODULE_PATH)
        self.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.module)

    def test_corrupt_or_truncated_part_is_rejected(self):
        for payload, size in [(b'abd', 3), (b'ab', 3)]:
            with self.subTest(payload=payload), self.assertRaises(ValueError):
                self.module.append_verified_part(io.BytesIO(payload), io.BytesIO(),
                    'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad',
                    size, hashlib.sha256())

    def test_parts_are_streamed_in_the_supplied_verified_order(self):
        output, combined = io.BytesIO(), hashlib.sha256()
        for block in (b'abc', b'def'):
            self.module.append_verified_part(io.BytesIO(block), output,
                                             hashlib.sha256(block).hexdigest(), 3, combined)
        self.assertEqual(output.getvalue(), b'abcdef')
        self.assertEqual(combined.hexdigest(), 'bef57ec7f53a6d40beb640a780a639c83bc29ac8a9816f1fc6c5c6dcd93c4721')

    def test_zip_traversal_and_symlink_are_rejected(self):
        for name in ('../panel', '/tmp/panel', 'a/../../panel', 'a\\panel'):
            with self.subTest(name=name), self.assertRaises(ValueError):
                self.module.safe_zip_member(zipfile.ZipInfo(name))
        symlink = zipfile.ZipInfo('work/panel')
        symlink.external_attr = 0o120777 << 16
        with self.assertRaises(ValueError):
            self.module.safe_zip_member(symlink)

    def test_extracted_member_must_match_its_own_source_digest(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / 'restored.parquet'
            archive = io.BytesIO()
            with zipfile.ZipFile(archive, 'w') as out:
                out.writestr('work/panel.parquet', b'abc')
            archive.seek(0)
            with zipfile.ZipFile(archive) as zin:
                member = zin.getinfo('work/panel.parquet')
                result = self.module.extract_verified_member(zin, member, target,
                    {'bytes': 3, 'sha256': 'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad'})
                self.assertEqual(target.read_bytes(), b'abc')
                self.assertEqual(result['bytes'], 3)
                with self.assertRaises(ValueError):
                    self.module.extract_verified_member(zin, member, target,
                        {'bytes': 3, 'sha256': '0' * 64})

    def test_duplicate_or_ambiguous_source_paths_are_rejected(self):
        records = [{'path': 'work/panel', 'bytes': 3, 'sha256': 'a'*64},
                   {'path': 'backup/work/panel', 'bytes': 3, 'sha256': 'b'*64}]
        with self.assertRaises(ValueError):
            self.module.find_manifest_record(records, 'work/panel')
        with self.assertRaises(ValueError):
            self.module.find_manifest_record(records, 'work/missing')

    def test_version_is_read_statically_without_unpickling(self):
        encoded = b'\x80\x04\x8c\x10_sklearn_version\x94\x8c\x051.7.2\x94.'
        self.assertEqual(self.module.static_sklearn_version(zlib.compress(encoded)), '1.7.2')
        self.assertIsNone(self.module.static_sklearn_version(zlib.compress(b'unrelated 1.2.3')))
        conflicting = encoded + b'\x8c\x10_sklearn_version\x94\x8c\x051.8.0\x94.'
        with self.assertRaises(ValueError):
            self.module.static_sklearn_version(zlib.compress(conflicting))

    def test_library_value_error_never_exports_its_private_message(self):
        self.assertTrue(callable(getattr(self.module, 'cli', None)),
                        'CLI must isolate approved diagnostics from library exception messages')
        output = io.StringIO()
        with (patch.object(self.module, 'main', side_effect=ValueError('PRIVATE_SOURCE_VALUE')),
              contextlib.redirect_stderr(output)):
            code = self.module.cli()
        self.assertEqual(code, 1)
        self.assertNotIn('PRIVATE_SOURCE_VALUE', output.getvalue())
        self.assertIn('ValueError', output.getvalue())


if __name__ == '__main__':
    unittest.main()
