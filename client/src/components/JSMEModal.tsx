import * as React from "react";
import { Modal } from "react-bootstrap";
import { Button } from "react-bootstrap";
// @ts-ignore
import { Jsme } from "jsme-react";

export function JSMEModal({
  open,
  initialSmiles = "",
  setOpen,
  onSave,
}: {
  open: boolean;
  initialSmiles?: string;
  setOpen?(show: boolean): void;
  onSave?(value: string | null): void;
}) {
  const [smiles, setSmiles] = React.useState<string>(initialSmiles);
  const [showJSME, setShowJSME] = React.useState<boolean>(false);
  const bodyRef = React.useRef<HTMLDivElement>(null);

  const handleClose = React.useCallback(() => {
    setOpen?.(false);
  }, [setOpen]);

  const handleSave = React.useCallback(() => {
    onSave?.(smiles || null);
  }, [onSave, smiles]);

  React.useEffect(() => {
    if (open) {
      const timeout = setTimeout(() => {
        setShowJSME(true);
      }, 200);
      return () => clearTimeout(timeout);
    } else {
      setShowJSME(false);
    }
  }, [open]);

  return (
    <>
      <Modal show={open} onHide={handleClose} size="lg">
        <Modal.Header closeButton>
          <Modal.Title>Draw molecule</Modal.Title>
        </Modal.Header>
        <Modal.Body>
          <div ref={bodyRef}></div>
          {showJSME ? (
            <Jsme
              height="600px"
              width={`${bodyRef.current?.getBoundingClientRect().width || 600}px`}
              options="star"
              // src="/jsme/jsme.nocache.js"
              // Load the distribution from the official "CDN".
              src="https://jsme-editor.github.io/dist/jsme/jsme.nocache.js"
              smiles={initialSmiles || undefined}
              onChange={setSmiles}
            />
          ) : null}
        </Modal.Body>
        <Modal.Footer>
          <Button variant="secondary" onClick={handleClose}>
            Close
          </Button>
          <Button variant="primary" onClick={handleSave}>
            Save
          </Button>
        </Modal.Footer>
      </Modal>
    </>
  );
}
