import * as React from "react";
import { Modal } from "react-bootstrap";
import { Button } from "react-bootstrap";

export function FileUploadModal({
  open,
  setOpen,
  onSave
}: {
  open: boolean;
  setOpen?(show: boolean): void;
  onSave?(value: string | null): void;
}) {
  const [input, setInput] = React.useState<string>('');

  const handleClose = React.useCallback(() => {
    setOpen?.(false);
  }, [setOpen]);

  const handleSave = React.useCallback(() => {
    onSave?.(input || null);
  }, [onSave, input]);

  return (
    <>
      <Modal show={open} onHide={handleClose} size="lg">
        <Modal.Header closeButton>
          <Modal.Title>Upload dataset</Modal.Title>
        </Modal.Header>
        <Modal.Body>
          <form>
          <div className="form-group">
            <label htmlFor="customFileInput">Either upload a file</label>
            <div className="custom-file">
              <input type="file" className="custom-file-input" id="customFileInput" onChange={(e) => {
                const files = Array.from(e.currentTarget.files || []);
                if(files.length === 1) {
                  files[0].text().then((res) => {
                    setInput(res);
                  })
                }
              }} />
              <label className="custom-file-label" htmlFor="customFileInput">Choose file</label>
            </div>
          </div>
          <div className="form-group">
            <label htmlFor="customFileTextarea">Or paste the content</label>
            <textarea id="customFileTextarea" className="form-control" onChange={(e) => setInput(e.currentTarget.value)} value={input.substr(0, 1000)} />
          </div>
          </form>
        </Modal.Body>
        <Modal.Footer>
          <Button variant="secondary" onClick={handleClose}>
            Close
          </Button>
          <Button variant="primary" onClick={handleSave}>
            Upload
          </Button>
        </Modal.Footer>
      </Modal>
    </>
  );
}
